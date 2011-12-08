package phylolab.mrp;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.Vector;
import phylolab.NewickTokenizer;

public class FastMRP {

	private static final String NEXUS = "NEXUS";
	private static final String NEXUS_HEADER = 	"#NEXUS\n" +
												"begin data;\n" +
													"\t dimensions ntax = @ nchar = @;\n" +
													"\t format missing = ?;\n" +
													"\tmatrix\n";
	private static final String NEXUS_FOOTER = "\t;\n" +
												"end;\n";
	private static final String PHYLIP = "PHYLIP";
	private static final String PHYLIP_HEADER = "@ @\n";
											

	char ONE = '1';
	char ZERO = '0';
	char MISSING = '?';	
	
	InfoPerTaxa perTaxaInfo = new InfoPerTaxa();
	TreeEndIndix treeEndIndix = new TreeEndIndix();
	
	// First list is the stack position, Each Collection is a bipartition (i.e index of taxa in one of the partitions)
	LinkedList<Collection<Integer>> stack = new LinkedList<Collection<Integer>>();
	
	int treeCount;
	
	private String treesFileName;
	private String mrpFileName;
	private String format;
	
	Random random = null;
	
	public FastMRP(String inFileName, String outFileName, String format) {
		treesFileName = inFileName;
		mrpFileName = outFileName;
		this.format = format;
	}
	
	public void setCharacters(char newOne, char newZero, char newMissing) {
		ONE = newOne;
		ZERO = newZero;
		MISSING = newMissing;
	}
	
	private void readTreesFile() throws IOException {
		BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(treesFileName)));		
		String tree;
		int treeId = 0;
		int lastColumnInd = -1;
	    //Read File Line By Line
	    while ((tree = in.readLine()) != null)   {
	    	NewickTokenizer tokenizer = new NewickTokenizer(tree);
	    	if (!"(".equals(tokenizer.nextToken())){
	    		throw new RuntimeException("The tree does not start with a (");
	    	}
	    	while (tokenizer.hasNext()) {
	    		String token = tokenizer.nextToken();
				if ("(".equals(token)){
					stack.addLast(new Vector<Integer>());
				} else if (")".equals(token)) {
					if (stack.size() > 0) {					
						Collection<Integer> top = stack.getLast();
						lastColumnInd = perTaxaInfo.addBipartition(top);
						stack.removeLast();
						if (stack.size() > 0) {
							stack.getLast().addAll(top);
						}
					}
				} else if (";".equals(token)) {
					treeEndIndix.addEndIndex(lastColumnInd);
				} else {
					Integer seqId = perTaxaInfo.mapNamesToIds(token);
					if (stack.size() > 0) {
						stack.getLast().add(seqId);
					}
					perTaxaInfo.addTreeToSequence(seqId, treeId);
				}
			}
	    	treeId++;
	    }
	    treeCount = treeId;
	}
	
	private void writeMRPToFile () throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(mrpFileName));
		writeHeaderToFile(out);
		// iterate over sequences, one by one
		perTaxaInfo.startSeqIteration();
		while (perTaxaInfo.hasMoreTaxa()) {
			perTaxaInfo.nextSequence();
			writeSequenceNameToFile(out);
			int column = 0;			
			Integer nextOneColumn = perTaxaInfo.nextBipartitionIndexForCurrentSequence(); 
			treeEndIndix.restartTreeIteration();
			for (int tree = 0; tree < treeCount; tree++) {
				int treeEndInd = treeEndIndix.getNexTreeEndInd();
				if (perTaxaInfo.currentSeqIsInTree(tree)) {
					while (column <= treeEndInd) {
						Boolean coding = perTaxaInfo.columnCoding.get(column);
						if (nextOneColumn == column) {
							out.write(coding? ONE:ZERO);
							nextOneColumn = perTaxaInfo.nextBipartitionIndexForCurrentSequence();
						} else {
							out.write(coding? ZERO:ONE);
						}
						column ++;
					}
				} else {
					int missingCount = treeEndInd - column + 1;
					char [] missings = new char[missingCount]; 
					Arrays.fill(missings, MISSING);
					out.write(missings);
					column = treeEndInd + 1;
				}
			}
			out.newLine();
		}
		writeFooterToFile(out);
		out.flush();
		out.close();
	}

	private void writeFooterToFile(BufferedWriter out) throws IOException {
		if (NEXUS.equalsIgnoreCase(format)){
			out.write(NEXUS_FOOTER);
		}		
	}

	private void writeHeaderToFile(BufferedWriter out) throws IOException {
		if (NEXUS.equalsIgnoreCase(format)){
			String header = NEXUS_HEADER.replaceFirst("@",perTaxaInfo.getNumberOfTaxa()+"").
				replaceFirst("@", perTaxaInfo.getNumberOfBipartitions()+"");
			out.write(header);
		} else if (PHYLIP.equalsIgnoreCase(format)){
			String header = PHYLIP_HEADER.replaceFirst("@",perTaxaInfo.getNumberOfTaxa()+"").
				replaceFirst("@", perTaxaInfo.getNumberOfBipartitions()+"");
			out.write(header);
		}
	}

	private void writeSequenceNameToFile(BufferedWriter out) throws IOException {
		if (NEXUS.equals(format)){
			out.write("\t'"+perTaxaInfo.getCurrentSeqName()+"' ");
		} else if (format.equalsIgnoreCase(PHYLIP)){
			out.write(perTaxaInfo.getCurrentSeqName()+" ");
		} else {
			out.write(">"+perTaxaInfo.getCurrentSeqName());
			out.newLine();
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 3) {
			System.err.println("Usage: <treesfile> <output> <ouputformat> [-dna] [-randomize seed]\n" +
					"		<treesfile>: A file containing newick trees, one tree per line\n" +
					"		<Output>: The name of the output MRP Matrix file\n" +
					"		<outformat>: use NEXUS for nexus, PHYLIP for phylip, or FASTA for fasta fromatted otuput\n"+
					" 		-dna: output As and Ts instead of 0 and 1\n" +
					"		-randomize: randomize 0-1 codings. Seed number is optional.");
					
			System.exit(1);
		}
		FastMRP mrpCon = new FastMRP(args[0], args[1],args[2]);
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-dna")) {
				mrpCon.setCharacters('A','T','-');
			} else if (args[i].equals("-randomize")) {
				try {
					int seed = Integer.parseInt(args[i+1]);
					mrpCon.random = new Random(seed);
				} catch(NumberFormatException e) {
					mrpCon.random = new Random();
				} catch(IndexOutOfBoundsException e) {
					mrpCon.random = new Random();
				}
			}
		}
		try {
			mrpCon.readTreesFile();
			mrpCon.writeMRPToFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	class InfoPerTaxa {
		Integer currentSeqIndex = -1;
		Integer currentColumnInd = -1;
		/*
		 * Index of 1st List ~ sequences, 2nd List's values ~ column indecies
		 */
		ArrayList<ArrayList<Integer>> columnsPerSequence = new ArrayList<ArrayList<Integer>>();
		/*
		 * Index of Arraylist ~ sequences, Elements in HashSet ~ Trees
		 */
		// TODO: can make this more memory efficient by turning the set to a list
		ArrayList<HashSet<Integer>> treesPerSequence = new ArrayList<HashSet<Integer>>();
		/*
		 * Mapping taxa names to IDs 
		 */
		HashMap<String, Integer> taxaNameToIndex = new HashMap<String, Integer>();
		ArrayList<String> taxaNames = new ArrayList<String>();
		/*
		 * Mapping columns to either 0 or 1 (randomly) to ensure equal 0s or 1s
		 * Index ~ 
		 */
		private ArrayList<Boolean> columnCoding = new ArrayList<Boolean>();
		
		private Iterator<String> namesIter;
		private Iterator<HashSet<Integer>> treesIter;
		private Iterator<ArrayList<Integer>> columnsIter;
		
		private Iterator<Integer> currentSeqColumnsIter;		
		private HashSet<Integer> currentTrees;		
		private ArrayList<Integer> currentColumns;
		private String currentSeqName;
						
		public Integer mapNamesToIds(String taxon) {
			if (taxaNameToIndex.containsKey(taxon)) {
				return taxaNameToIndex.get(taxon);
			} else {
				currentSeqIndex++;
				taxaNameToIndex.put(taxon, currentSeqIndex);
				taxaNames.add(currentSeqIndex,taxon);
				return currentSeqIndex;
			}
		}
		
		public int addBipartition(Collection<Integer> bipartition) {
			currentColumnInd++;
			// Randomly assign codings to columns if -randomize is provided
			columnCoding.add(currentColumnInd , random==null? true:random.nextBoolean()); 

			for (Integer sequence : bipartition) {
				columnsPerSequence.get(sequence).add(currentColumnInd);
			}
			return currentColumnInd;
		}

		void addTreeToSequence(Integer seq, Integer tree){
			// If the sequence is encountered for the first time, add it to datastructures.
			if (seq >= treesPerSequence.size()) {
				treesPerSequence.add(seq,new HashSet<Integer>());
				columnsPerSequence.add(seq, new ArrayList<Integer>());
			}
			treesPerSequence.get(seq).add(tree);
		}
		
		void startSeqIteration() {
			treesIter = treesPerSequence.iterator();
			columnsIter = columnsPerSequence.iterator();
			namesIter = taxaNames.iterator();
		}
		
		boolean hasMoreTaxa () {
			return treesIter.hasNext();
		}
		
		void nextSequence () {
			currentTrees = treesIter.next();
			currentColumns = columnsIter.next();
			currentSeqColumnsIter = currentColumns.iterator();
			currentSeqName = namesIter.next();
		}
		boolean currentSeqIsInTree (Integer treeId) {
			return currentTrees.contains(treeId);
		}
		public Integer nextBipartitionIndexForCurrentSequence() {
			if (currentSeqColumnsIter.hasNext()) {
				return currentSeqColumnsIter.next();
			} else {
				return -1;
			}
		}
		public String getCurrentSeqName() {
			return currentSeqName;
		}
		
		public Integer getNumberOfTaxa(){
			return currentSeqIndex+1;
		}
		public Integer getNumberOfBipartitions(){
			return currentColumnInd+1;
		}
	}	
	
	class TreeEndIndix {
		/*
		 * Index to list~ tree, Integer ~ the index of the last column of the tree
		 */
		ArrayList<Integer> treeEndIndex = new ArrayList<Integer>();
		private Iterator<Integer> iter;
		public void addEndIndex(int lastColumnIndex) {
			treeEndIndex.add(lastColumnIndex);
		}
		
		void restartTreeIteration(){
			iter = treeEndIndex.iterator();
		}
		
		int getNexTreeEndInd(){
			return iter.next();
		}
	}
}
