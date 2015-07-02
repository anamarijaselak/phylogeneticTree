package hr.fer.phylogeny;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class PhylogeneticTree {

	static char[] nucleotides = new char[4];
	static int brTopology = 0;

	static List<Integer> informativePositions = new ArrayList<Integer>();
	static List<FullTree> fullTrees = new ArrayList<FullTree>();
	static int brPos;
	static int addCost = 0;
	static int minimumCost = 0;
	static int brStabala = 0;
	static int pom = 0;
	static Map<Integer, String> sequences = new HashMap<Integer, String>();

	public static void main(String[] args) {

		if (args.length != 2) {
			System.out.println("<input file> <output directory>");
			System.exit(0);
		}

		Path fastaPath = Paths.get(args[0]);
		String directory = args[1];
		List<String> sequences = new ArrayList<String>();
		List<String> lines = new ArrayList<String>();

		try {
			lines = Files.readAllLines(fastaPath);
		} catch (IOException e) {
			System.out.println("Unable to read config file!");

		}
		parseClustalWFile(lines, sequences);
		if (sequences.size() <= 3) {
			System.out.println("File must contain at least 4 sequences");
			System.exit(0);
		}

		nucleotides[0] = 'A';
		nucleotides[1] = 'T';
		nucleotides[2] = 'C';
		nucleotides[3] = 'G';

		List<TreeNode> trees = null;

		brPos = sequences.get(0).length();

		getInformativePositions(sequences, brPos);

		Map<Integer, List<TreeNode>> treesPerPosition = new HashMap<Integer, List<TreeNode>>();
		Map<Integer, List<TopologyTrees>> almostFullTreesPerPos = new HashMap<Integer, List<TopologyTrees>>();
		Set<Integer> topologies = new LinkedHashSet<Integer>();

		brTopology = 0;
		int brTop = 1;
		int seqSize;
		if (sequences.size() == 4) {
			seqSize = 5;
		} else {
			seqSize = sequences.size();
		}

		System.out.println("Constructing optimal tree..");
		for (int brSeq = 3; brSeq < seqSize - 1; brSeq++) {

			int brConf = 2 * (brSeq) - 3;
			brTop = brTop * brConf;

			for (Integer position : informativePositions) {

				if (brSeq == 3) {
					trees = new ArrayList<TreeNode>();
					List<TreeNode> unrootedTrees = new ArrayList<TreeNode>();
					trees = new ArrayList<TreeNode>();
					Node n1 = new Node(sequences.get(0).charAt(position), null,
							null, null);
					n1.brSeq = 0;
					Node n2 = new Node(sequences.get(1).charAt(position), null,
							null, null);
					n2.brSeq = 1;
					Node n3 = new Node(sequences.get(2).charAt(position), null,
							null, null);
					n3.brSeq = 2;

					for (int i = 0; i < 4; i++) {
						unrootedTrees.add(new Node(nucleotides[i], n1, n2, n3));
					}

					// dodavanje cetvrtog i rooting pocinje
					PhylogeneticTree.rooting(unrootedTrees, trees, position,
							sequences.get(3));
					treesPerPosition.put(position, trees);

					if (sequences.size() == 4) {

						continue;
					}

					if (brSeq == sequences.size() - 2) {

						almostFullTreesPerPos.put(position,
								getTreesForTopologies(trees, topologies));
					}
					continue;

				}

				Topology top = new Topology(brTopology);

				List<TreeNode> newTrees = new ArrayList<TreeNode>();
				Node addNode = new Node(sequences.get(brSeq).charAt(position),
						null, null, null);
				addNode.brSeq = brSeq;

				int size = treesPerPosition.get(position).size();
				// System.out.println("Velicina "+size+" brseq "+(brSeq+1));
				for (int brTree = 0; brTree < size; brTree++) {

					TreeNode root = treesPerPosition.get(position).get(brTree);
					top.brTopology = brTopology
							+ (brTopology - ((Root) root).topology) * brConf;

					makeNewTrees(root, addNode, (Root) root, newTrees, top);

				}
				treesPerPosition.put(position, newTrees);

				if (brSeq == sequences.size() - 2) {

					almostFullTreesPerPos.put(position,
							getTreesForTopologies(newTrees, topologies));
				}

			}
			brTopology += brTop;
			if (sequences.size() == 4) {
				calculateForFour(treesPerPosition);
			}

		}

		if (sequences.size() > 4) {
			int i = 0;
			int brConf = 2 * (sequences.size() - 1) - 3;
			Topology top = new Topology(brTopology);

			for (Integer topology : topologies) {

				Map<Integer, List<TreeNode>> lastTrees = new HashMap<Integer, List<TreeNode>>();
				i++;
				if (i == 1) {
					for (Integer j : informativePositions) {

						Node addNode = new Node(sequences.get(
								sequences.size() - 1).charAt(j), null, null,
								null);
						addNode.brSeq = sequences.size() - 1;

						buildFullPos(addNode, almostFullTreesPerPos, j, top,
								lastTrees, topology, brConf);

					}
					pom++;
					buildFullTrees(lastTrees);

					minimumCost = getBoundary();

					continue;

				}

				// dobili smo prvu granicu, sada prije nego dodamo zadnju
				// sekvencu
				// moramo izracunati cijenu ukupnog stabla u prethodnom koraku
				// za
				// određenu topologiju
				// ako je cijena već veća tu topologiju preskačemo

				int almostFullCost = 0;
				for (Integer j : informativePositions) {

					for (TopologyTrees optTop : almostFullTreesPerPos.get(j)) {
						if (optTop.topology == topology) {
							almostFullCost += optTop.getMinimalCost();
						}
					}

				}

				if (almostFullCost > minimumCost) {

					continue;
				}

				// cijena prethodne razinne je manja ili jednaka dosadasnjoj
				// granici
				// pa za ovu topologiju gradimo sva stabla
				// izgrađena stabla sada ostavljaju samo ona optimalna za svaku
				// topologiju
				// spajaju se i spremaju u fullTrees

				for (Integer j : informativePositions) {

					Node addNode = new Node(sequences.get(sequences.size() - 1)
							.charAt(j), null, null, null);
					addNode.brSeq = sequences.size() - 1;

					buildFullPos(addNode, almostFullTreesPerPos, j, top,
							lastTrees, topology, brConf);

				}

				pom++;
				buildFullTrees(lastTrees);

				int cost = getBoundary();

				if (cost < minimumCost) {
					minimumCost = cost;

				}

			}

		}

		System.out.println("Optimal tree cost : " + minimumCost);
		int brT = 0;
		FullTree[] optimalTrees = new FullTree[fullTrees.size()];

		for (FullTree fullTree : fullTrees) {

			if (fullTree.totalCost == minimumCost) {

				optimalTrees[brT] = fullTree;
				brT++;

				String fileName = directory + "optimalTree" + brT + ".nwk";
				Path path = Paths.get(fileName);
				dumpToGraph(fullTree.representative, path);

			}
		}
		System.out.println("Number of optimal tree topologies : " + brT);

	}

	public static void calculateForFour(
			Map<Integer, List<TreeNode>> treesPerPosition) {
		buildFullTrees(treesPerPosition);
		minimumCost = getBoundary();

	}

	public static void parseClustalWFile(List<String> lines,
			List<String> sequences2) {
		Map<String, String> sequences = new LinkedHashMap<String, String>();
		lines.remove(0);
		boolean started = false;

		for (String line : lines) {
			if (line.equals("") && !started) {
				continue;
			}
			if (!line.equals("") && !started) {
				started = true;
			}
			if (line != null && !line.equals("") && !line.isEmpty()) {
				if (!line.matches(".*[a-zA-Z]+.*")) {

					continue;
				}
				String[] elements = line.replaceAll("\\s+", " ").split(" ");

				if (elements.length == 0) {
					continue;
				}
				if (sequences.containsKey(elements[0])) {
					String element = sequences.get(elements[0]).concat(
							elements[1]);
					sequences.put(elements[0], element);

				} else {
					sequences.put(elements[0], elements[1]);

				}

			}

		}
		int br = 0;

		for (Map.Entry<String, String> entry : sequences.entrySet()) {

			// System.out.println(entry.getKey());
			PhylogeneticTree.sequences.put(br, entry.getKey());

			sequences2.add(entry.getValue());

			br++;

		}

	}

	/**
	 * Metoda koja izgraduje cijelo stablo, spaja sve pozicije. Stabla koja
	 * prima su stabla za sve pozicije za određenu topologiju.
	 * 
	 * @param lastTrees
	 * @param fullTrees
	 */
	private static void buildFullTrees(Map<Integer, List<TreeNode>> lastTrees) {

		Set<Integer> topologies = new HashSet<Integer>();
		// dodajem sve topologije

		for (TreeNode tree : lastTrees.get(informativePositions.get(0))) {

			topologies.add(((Root) tree).topology);

		}

		// mapa u kojoj imamo kao kljuc poziciju, a kao vrijednost listu svih
		// klasa koje predstavljaju topologiju + stabla koja pripadaju toj
		// topologiji
		// U mapi ce se nalaziti samo ona stabla koja su optimalna za neku
		// topologiju jer se zove metoda optimize
		Map<Integer, List<TopologyTrees>> optTreeMap = new HashMap<Integer, List<TopologyTrees>>();

		for (Integer i : informativePositions) {

			List<TopologyTrees> treePos = new ArrayList<TopologyTrees>();

			for (TreeNode tree : lastTrees.get(i)) {

				boolean found = false;
				for (TopologyTrees optTree : treePos) {
					if (((Root) tree).topology == optTree.topology) {
						optTree.addTree(tree);
						found = true;
					}
				}
				if (!found) {
					TopologyTrees optTree = new TopologyTrees(
							((Root) tree).topology, new ArrayList<TreeNode>());
					optTree.position = i;
					optTree.addTree(tree);

					treePos.add(optTree);
				}
			}

			for (TopologyTrees optTop : treePos) {

				// System.out.println("Zovem optimize za "+i);
				optTop.optimize();
			}

			optTreeMap.put(i, treePos);
		}

		// izgradimo sva stabla
		// u fullTrees imamo stabla za sve konfiguracije nastale od stabla s
		// jednom sekv manje ali od tocno odredene topologije
		for (Integer topology : topologies) {

			FullTree fullTree = new FullTree();
			TreeNode representative = null;
			int totalCost = 0;

			for (Integer i : informativePositions) {

				// idem za svaku poziciju i uzimam cijenu koja je vec zapisana
				// jer sam prije pozvala metodu optimize
				boolean haveRepresentative = false;
				for (TopologyTrees optTop : optTreeMap.get(i)) {

					if (optTop.topology == topology) {

						totalCost += optTop.cost;
						// System.out.println("Pozicija "+i+" cijena "+optTop.cost);

						if (!haveRepresentative) {

							representative = optTop.trees.get(0);
							haveRepresentative = true;
						}
					}
				}
			}

			fullTree.representative = representative;

			fullTree.topology = topology;
			fullTree.totalCost = totalCost + addCost;

			fullTrees.add(fullTree);

		}

	}

	/**
	 * Za određenu topologiju dodaje zadnju sekvencu. Vraća listu stabala za tu
	 * poziciju i tu topologiju.
	 * 
	 * @param addNode
	 * @param optimalTressPos
	 * @param j
	 * @param top
	 * @param lastTrees
	 * @param topology
	 * @param brConf
	 */
	private static void buildFullPos(Node addNode,
			Map<Integer, List<TopologyTrees>> optimalTressPos, int j,
			Topology top, Map<Integer, List<TreeNode>> lastTrees,
			Integer topology, int brConf) {
		List<TreeNode> newTrees = new ArrayList<TreeNode>();
		List<TreeNode> topTrees = null;

		for (TopologyTrees optTop : optimalTressPos.get(j)) {
			if (optTop.topology == topology) {
				topTrees = optTop.trees;
			}
		}

		for (TreeNode root : topTrees) {

			top.brTopology = brTopology + (brTopology - topology) * brConf;
			makeNewTrees(root, addNode, (Root) root, newTrees, top);

		}
		lastTrees.put(j, newTrees);

	}

	private static int getBoundary() {

		int minimum = fullTrees.get(0).totalCost;
		for (FullTree fullTree : fullTrees) {
			if (fullTree.totalCost < minimum) {
				minimum = fullTree.totalCost;
			}
		}
		return minimum;

	}

	/**
	 * Metoda vraća listu objekata tipa Topology. Razvrstava sva stabla s
	 * obzirom na topologije.
	 * 
	 * @param trees
	 * @param topologies
	 * @return
	 */
	private static List<TopologyTrees> getTreesForTopologies(
			List<TreeNode> trees, Set<Integer> topologies) {
		List<TopologyTrees> treePos = new ArrayList<TopologyTrees>();
		for (TreeNode tree : trees) {
			topologies.add(((Root) tree).topology);
			boolean found = false;
			for (TopologyTrees optTree : treePos) {
				if (((Root) tree).topology == optTree.topology) {
					optTree.addTree(tree);
					found = true;
				}
			}
			if (!found) {
				TopologyTrees optTree = new TopologyTrees(
						((Root) tree).topology, new ArrayList<TreeNode>());
				optTree.addTree(tree);
				treePos.add(optTree);
			}
		}

		return treePos;
	}

	private static void getInformativePositions(List<String> sequences,
			int seqLength) {

		for (int i = 0; i < seqLength; i++) {
			int[] brN = new int[4];
			boolean flag_add = true;
			for (String s : sequences) {

				if (s.charAt(i) == 'A') {
					brN[0]++;

					if (brN[1] != 0 || brN[2] != 0 || brN[3] != 0) {
						flag_add = false;

					}

				} else if (s.charAt(i) == 'T') {
					brN[1]++;

					if (brN[0] != 0 || brN[2] != 0 || brN[3] != 0) {
						flag_add = false;

					}

				} else if (s.charAt(i) == 'C') {
					brN[2]++;

					if (brN[1] != 0 || brN[0] != 0 || brN[3] != 0) {
						flag_add = false;

					}

				} else if (s.charAt(i) == 'G') {
					brN[3]++;

					if (brN[1] != 0 || brN[2] != 0 || brN[0] != 0) {
						flag_add = false;

					}

				}
			}
			if (brN[0] == sequences.size() - 1
					|| brN[1] == sequences.size() - 1
					|| brN[2] == sequences.size() - 1
					|| brN[3] == sequences.size() - 1) {

				flag_add = true;
				if (brN[0] == 1 || brN[1] == 1 || brN[2] == 1 || brN[3] == 1) {
					addCost++;
				}
			}
			if (!flag_add) {

				informativePositions.add(i);
			}

		}

	}

	private static void rooting(List<TreeNode> unrootedTrees,
			List<TreeNode> trees, int position, String seq) {
		for (int i = 0; i < unrootedTrees.size(); i++) {
			Node n = (Node) unrootedTrees.get(i);
			for (int j = 0; j < 4; j++) {
				Node left, right, parent;
				Node n4 = new Node(seq.charAt(position), null, null, null);
				n4.brSeq = 3;
				left = n.leftChild.copy();
				right = n4.copy();
				parent = n.copy();
				Node newNode1 = new Node(nucleotides[j], left, n4.copy(),
						parent);

				left.parent = newNode1;
				right.parent = newNode1;
				parent.leftChild = n.rightChild.copy();
				parent.rightChild = n.parent.copy();
				parent.parent = newNode1;

				Root r = new Root(newNode1, parent, 1);

				trees.add(r);

				left = n.rightChild.copy();
				right = n4.copy();
				parent = n.copy();
				Node newNode2 = new Node(nucleotides[j], right, left, parent);

				left.parent = newNode2;
				right.parent = newNode2;

				parent.leftChild = n.leftChild.copy();
				parent.rightChild = n.parent.copy();
				parent.parent = newNode2;

				trees.add(new Root(newNode2, parent, 2));

				left = n.parent.copy();
				right = n4.copy();
				parent = n.copy();
				Node newNode3 = new Node(nucleotides[j], right, left, parent);

				left.parent = newNode3;
				right.parent = newNode3;

				parent.leftChild = n.leftChild.copy();
				parent.rightChild = n.rightChild.copy();
				parent.parent = newNode3;

				trees.add(new Root(newNode3, parent, 3));

			}

		}

	}

	public static void dumpToGraph(TreeNode root, Path path) {

		String newickFormat = inOrderNewick(root);
		List<String> lines = new ArrayList<String>();
		lines.add(newickFormat);
		try {

			Files.write(path, lines);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static String inOrderNewick(TreeNode root) {
		if (root.hasChild()) {
			String output = "";
			output += "(";
			output += inOrderNewick(root.leftChild);
			output += ",";
			output += inOrderNewick(root.rightChild);
			output += ")";
			return output;
		} else {
			return String.valueOf(sequences.get(root.brSeq));
		}
	}

	public static void printBinaryTree(TreeNode root, int level) {
		if (root == null)
			return;
		printBinaryTree(root.rightChild, level + 1);
		if (level != 0) {
			for (int i = 0; i < level - 1; i++)
				System.out.print("|\t");
			System.out.println("|-------" + sequences.get(root.brSeq));
		} else
			System.out.println(sequences.get(root.brSeq));
		printBinaryTree(root.leftChild, level + 1);

	}

	public static void makeNewTrees(TreeNode currentNode, Node addNode,
			Root root, List<TreeNode> roots, Topology top) {

		if (currentNode == null) {
			return;
		}

		// gradim novo stablo tako da stvaram novi node, dodajem ga kao lijevo
		// dijete na trenutni cvor
		// provjeri da ne dodajes na listove
		if (currentNode.leftChild != null && currentNode.rightChild != null) {
			top.brTopology++;
			// dodajem zmeđu roditelja i lijevog djeteta
			for (int i = 0; i < 4; i++) {
				Node left = currentNode.leftChild;
				Node right = addNode;
				Node parent;
				if (currentNode instanceof Root) {
					parent = currentNode.rightChild;

				} else {
					parent = ((Node) currentNode);
				}

				Node newNode = new Node(nucleotides[i], left, right, parent);
				if (currentNode instanceof Root) {
					parent.parent = newNode;
				}
				currentNode.leftChild = newNode;

				// pozovi za kopiranje -> posalji samo root
				Root newRoot = new Root();

				newRoot.topology = top.brTopology;
				newRoot.nucleotide = 'r';
				PhylogeneticTree.copyTree(root, newRoot);
				roots.add(newRoot);

				// vrati stablo na staro

				currentNode.leftChild = newNode.leftChild;
			}

			// dodajem izmedu roditelja i desnog djeteta osim ako roditelj nije
			// korijen

			if (!(currentNode instanceof Root)) {
				top.brTopology++;
				for (int i = 0; i < 4; i++) {
					Node left = currentNode.rightChild;
					Node right = addNode;
					Node parent;

					parent = ((Node) currentNode);

					Node newNode = new Node(nucleotides[i], left, right, parent);

					currentNode.rightChild = newNode;

					// pozovi za kopiranje -> posalji samo root
					Root newRoot = new Root();
					// System.out.println("Top -> " + top.brTopology);
					newRoot.topology = top.brTopology;
					newRoot.nucleotide = 'r';

					PhylogeneticTree.copyTree(root, newRoot);
					roots.add(newRoot);

					// vrati stablo na staro

					currentNode.rightChild = newNode.leftChild;
				}

			}

		}

		makeNewTrees(currentNode.leftChild, addNode, root, roots, top);

		makeNewTrees(currentNode.rightChild, addNode, root, roots, top);

	}

	public static void copyTree(TreeNode oldNode, TreeNode newNode) {

		if (oldNode.leftChild != null) {
			newNode.leftChild = oldNode.leftChild.copy();
			copyTree(oldNode.leftChild, newNode.leftChild);
		}

		if (oldNode.rightChild != null) {
			newNode.rightChild = oldNode.rightChild.copy();
			copyTree(oldNode.rightChild, newNode.rightChild);
		}

	}

	public static int countTreeCost(TreeNode currentNode) {

		if (currentNode == null) {
			return 0;
		}

		int cost = 0;
		if (currentNode instanceof Root) {
			cost += getCostBetween(currentNode.rightChild,
					currentNode.leftChild);

		} else {
			if (currentNode.leftChild != null) {
				cost += getCostBetween(currentNode, currentNode.leftChild);
			}
			if (currentNode.rightChild != null) {
				cost += getCostBetween(currentNode, currentNode.rightChild);
			}

		}

		return cost + countTreeCost(currentNode.leftChild)
				+ countTreeCost(currentNode.rightChild);

	}

	private static int getCostBetween(TreeNode first, TreeNode second) {
		if (first.nucleotide != second.nucleotide && first.nucleotide != '-'
				&& second.nucleotide != '-') {
			return 1;
		} else {
			return 0;
		}

	}

}
