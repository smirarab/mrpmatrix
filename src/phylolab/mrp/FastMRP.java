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

	final char ONE = '1';
	final char ZERO = '0';
	final char MISSING = '?';	
	
	InfoPerTaxa perTaxaInfo = new InfoPerTaxa();
	TreeEndIndix treeEndIndix = new TreeEndIndix();
	
	// First list is the stack position, Each Collection is a bipartition (i.e index of taxa in one of the partitions)
	LinkedList<Collection<Integer>> stack = new LinkedList<Collection<Integer>>();
	
	int treeCount;
	
	private String treesFileName;
	private String mrpFileName;
	private String format;
	
	public FastMRP(String inFileName, String outFileName, String format) {
		treesFileName = inFileName;
		mrpFileName = outFileName;
		this.format = format;
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
						if (nextOneColumn == column) {
							out.write(ONE);
							nextOneColumn = perTaxaInfo.nextBipartitionIndexForCurrentSequence();
						} else {
							out.write(ZERO);
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
		if (NEXUS.equals(format)){
			out.write(NEXUS_FOOTER);
		}		
	}

	private void writeHeaderToFile(BufferedWriter out) throws IOException {
		if (NEXUS.equals(format)){
			String header = NEXUS_HEADER.replaceFirst("@",perTaxaInfo.getNumberOfTaxa()+"").
				replaceFirst("@", perTaxaInfo.getNumberOfBipartitions()+"");
			out.write(header);
		}		
	}

	private void writeSequenceNameToFile(BufferedWriter out) throws IOException {
		if (NEXUS.equals(format)){
			out.write("\t'"+perTaxaInfo.getCurrentSeqName()+"' ");
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
			System.err.println("Usage: <treesfile> <output> <ouputformat>\n" +
					"		<treesfile>: A file containing newick trees, one tree per line\n" +
					"		<Output>: The name of the output MRP Matrix file\n" +
					"		<outformat>: use NEXUS for nexus, or FASTA for fasta fromatted otuput\n");
			System.exit(1);
		}
		FastMRP mrpCon = new FastMRP(args[0], args[1],args[2]);
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
		 * 1st List's index ~ sequences, 2nd List's values ~ column indecies
		 */
		ArrayList<ArrayList<Integer>> columnsPerSequence = new ArrayList<ArrayList<Integer>>();
		/*
		 * Index to Arraylist ~ sequences, Elements in HashSet ~ Trees
		 */
		ArrayList<HashSet<Integer>> treesPerSequence = new ArrayList<HashSet<Integer>>();
		/*
		 * Mapping taxa names to IDs 
		 */
		HashMap<String, Integer> taxaNameToIndex = new HashMap<String, Integer>();
		
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
				return currentSeqIndex;
			}
		}
		
		public int addBipartition(Collection<Integer> bipartition) {
			currentColumnInd++;
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
			namesIter = taxaNameToIndex.keySet().iterator();
		}		
		boolean hasMoreTaxa () {
			return treesIter.hasNext();
		}
		void nextSequence () {
			currentTrees = treesIter.next();
			currentColumns = columnsIter.next();
			currentSeqName = namesIter.next();
			currentSeqColumnsIter = currentColumns.iterator();			
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