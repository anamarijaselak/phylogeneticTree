package hr.fer.phylogeny;

public class Node extends TreeNode {

	public Node parent;

	public Node(char nucleotide, Node leftChild, Node rightChild, Node parent) {
		super(leftChild, rightChild, nucleotide);

		this.parent = parent;

	}

	public Node() {

	}

	public Node copy() {
		Node newN = new Node(nucleotide, leftChild, rightChild, parent);
		newN.brSeq = this.brSeq;
		return newN;
	}

}
