from rnadist.utils import list_internal_leafsets, Node

def test_list_internal_leafsets():

    L1 = Node(label=1) 
    L2 = Node(label=2) 
    L3 = Node(label=3) 
    L4 = Node(label=4)
    L5 = Node(label=5)
    L6 = Node(label=6)

    R = Node()
    I1 = Node()
    I2 = Node()

    R.add_child(L1)
    R.add_child(L6)
    R.add_child(I1)
    R.add_child(I2)

    I1.add_child(L2)
    I1.add_child(L3)
    
    I2.add_child(L4)
    I2.add_child(L5)
    
    assert(list_internal_leafsets(R)==[[1,6],[2,3],[4,5]])
