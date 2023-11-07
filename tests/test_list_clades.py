from rnadist.utils import Node, list_clades

def test_list_clades():
    root = Node()
    I1 = Node()
    I2 = Node()
    L1 = Node(label=1)
    L2 = Node(label=2)
    L3 = Node(label=3)
    L4 = Node(label=4)
    root.add_child(I1)
    root.add_child(I2)
    I1.add_child(L1)
    I1.add_child(L2)
    I2.add_child(L3)
    I2.add_child(L4)

    assert(list_clades(root)==[[1],[2],[1,2],[3],[4],[3,4],[1,2,3,4]])
