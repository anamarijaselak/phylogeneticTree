package hr.fer.phylogeny;

public class UnrootedTree {

	private Node internal;
	private Node child1;
	private Node child2;
	private Node child3;

	public UnrootedTree(Node internal, Node child1, Node child2, Node child3) {
		super();
		this.internal = internal;
		this.child1 = child1;
		this.child2 = child2;
		this.child3 = child3;
	}

	public Node getInternal() {
		return internal;
	}

	public void setInternal(Node internal) {
		this.internal = internal;
	}

	public Node getChild1() {
		return child1;
	}

	public void setChild1(Node child1) {
		this.child1 = child1;
	}

	public Node getChild2() {
		return child2;
	}

	public void setChild2(Node child2) {
		this.child2 = child2;
	}

	public Node getChild3() {
		return child3;
	}

	public void setChild3(Node child3) {
		this.child3 = child3;
	}

}
