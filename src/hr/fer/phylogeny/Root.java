package hr.fer.phylogeny;

public class Root extends TreeNode {

	public int topology;

	public Root(Node leftChild, Node rightChild, int topology) {
		super(leftChild, rightChild, 'r');

		this.topology = topology;

	}

	public Root() {

	}

}
