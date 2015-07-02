package hr.fer.phylogeny;

import java.util.ArrayList;
import java.util.List;

public class TopologyTrees {

	public int topology;
	public List<TreeNode> trees;
	public int position;
	public int cost;

	public TopologyTrees(int top, List<TreeNode> trees) {
		this.topology = top;
		this.trees = trees;
	}

	public void addTree(TreeNode tree) {
		this.trees.add(tree);
	}

	public void optimize() {

		List<TreeNode> forRemove = new ArrayList<TreeNode>();
		int minimumCost = PhylogeneticTree.countTreeCost(trees.get(0));

		for (TreeNode tree : trees) {
			int cost = PhylogeneticTree.countTreeCost(tree);

			if (cost < minimumCost) {
				minimumCost = cost;
			}
		}

		for (TreeNode tree : trees) {

			if (PhylogeneticTree.countTreeCost(tree) > minimumCost) {
				forRemove.add(tree);
			}
		}
		this.cost = minimumCost;
		// System.out.println("Cijena je "+this.cost);
		trees.removeAll(forRemove);

	}

	@Override
	public boolean equals(Object arg) {
		TopologyTrees other = (TopologyTrees) arg;
		if (this.topology == other.topology) {
			return true;
		}
		return false;
	}

	public int getMinimalCost() {

		int minimumCost = PhylogeneticTree.countTreeCost(trees.get(0));
		for (TreeNode tree : trees) {
			int cost = PhylogeneticTree.countTreeCost(tree);
			if (cost < minimumCost) {
				minimumCost = cost;
			}
		}
		return minimumCost;
	}

}
