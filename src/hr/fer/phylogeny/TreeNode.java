package hr.fer.phylogeny;

public abstract class TreeNode {

	public Node leftChild;
	public Node rightChild;
	public char nucleotide;
	public int brSeq;

	public TreeNode(Node leftChild, Node rightChild, char nucleotide) {
		super();
		this.leftChild = leftChild;
		this.rightChild = rightChild;
		this.nucleotide = nucleotide;
	}

	public TreeNode() {
		// TODO Auto-generated constructor stub
	}

	public boolean hasChild() {
		if (leftChild == null && rightChild == null) {
			return false;
		}
		return true;
	}

}
