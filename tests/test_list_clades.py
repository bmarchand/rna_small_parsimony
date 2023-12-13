from rnadist.utils import Node, list_clades, RFdistance

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

    L = [[1,2],[3,4],[1,2,3,4]]
    res = [frozenset(s) for s in L]

    print(list_clades(root))
    assert(list_clades(root)==res)

def test_rf_distance():

    root1 = Node()
    I1 = Node()
    I2 = Node()
    L1 = Node(label=1)
    L2 = Node(label=2)
    L3 = Node(label=3)
    L4 = Node(label=4)
    root1.add_child(I1)
    root1.add_child(I2)
    I1.add_child(L1)
    I1.add_child(L2)
    I2.add_child(L3)
    I2.add_child(L4)
    
    root2 = Node()
    I12 = Node()
    I22 = Node()
    L12 = Node(label=1)
    L22 = Node(label=3)
    L32 = Node(label=2)
    L42 = Node(label=4)
    root2.add_child(I12)
    root2.add_child(I22)
    I12.add_child(L12)
    I12.add_child(L22)
    I22.add_child(L32)
    I22.add_child(L42)

    assert(RFdistance(root1, root1)==0)
    assert(RFdistance(root2, root2)==0)
    assert(RFdistance(root1, root2)==4)
