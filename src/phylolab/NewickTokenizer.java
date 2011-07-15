package phylolab;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class NewickTokenizer {

	String tree;
	Matcher pattern;
	boolean strip;
	
	public NewickTokenizer (String input) {
		init(input, true);
	}
	
	public NewickTokenizer (String input, boolean strip) {
		init(input,strip);		
	}

	private void init(String input,boolean strip) {
		tree = input;
		this.strip = strip;
		if (strip) {
			pattern = Pattern.compile("([(])|([)][^,:;)]*)|([;])|(:)|([^,);(:]*)").matcher(tree);
		} else {
			pattern = Pattern.compile("([(])|([)][^,;)]*)|([;])|([^,);(]*)").matcher(tree);
		}
		pattern.find();
	}
	
	public boolean hasNext() {
		return !pattern.hitEnd();
	}
	public String nextToken() {
		String res = pattern.group();
		pattern.find();
		// This is to strip off any support value / internal label nodes that can follow a left bracket.
		if (strip && res.startsWith(")")) {
			return ")";
		}
		if ("".equals(res)) {
			return nextToken();
		} else if (":".equals(res)){
			// This is to strip off the branch lenght values.
			nextToken();
			return nextToken();
		} else {
			return res;
		}		
	}
	
	public static void main (String [] args) {
		//String tree = "(t90:0.01,((t18:0.04,(((t17:0.0,t49:0.09):0.00,(((((t81:0.037,t45:0.02):0.02,t66:0.01):0.00724936580827239413,t41:0.05336582750396093311):0.03045292104932446550,(t14:0.11614645744250107207,((t16:0.05973131223865356387,t91:0.06115882420120943158):0.03367781109364907655,((t30:0.06168112199284891961,(t82:0.00707911640722609301,t19:0.01702795908063739483):0.04105036384987210962):0.04367847525636570777,(t88:0.08520512585491028801,((t77:0.02231272403400484661,t100:0.03396176535082153641):0.03970879193218978392,t46:0.07474947196269035588):0.03491872410661160664):0.00983925300854607623):0.00000122920446704147):0.00271415429521279158):0.00799815476646543837):0.01024163042108295132,(t80:0.07155782866412271903,t68:0.09571882482752180898):0.01525096682579810646):0.01622216107791855585):0.00574275591113855306,t31:0.08260865518257545781):0.00767107143981279709):0.00100486256845222326,t15:0.08070670351238432016):0.03897657131187436119,t43:0.01734246018246490828):0.0;";
		String tree = "(2sd:32[2],((M:29[23],(A:10[2],(C:1,D:3)0.30[12]:232,B:12),N)22:232[12],L));";
		NewickTokenizer tokenizer = new NewickTokenizer(tree,false);
		while (tokenizer.hasNext()) {
			System.out.print(tokenizer.nextToken()+"  ");			
		}
	}
	
}

