from rnadist.utils import PhyloNode, fitch_with_trait_vector, VectorTraits 
import numpy as np

def test_small_tree1():

    Root = PhyloNode('root', False)
    I1 = PhyloNode('I1', False)
    L1 = PhyloNode('L1', True)
    L2 = PhyloNode('L2', True)
    L3 = PhyloNode('L3', True)

    Root.add_children(I1,L3)
    I1.add_children(L1,L2)

    vector_traits = VectorTraits(3)
    vector_traits.add_trait('L1', [0,1,0])
    vector_traits.add_trait('L2', [1,1,0])
    vector_traits.add_trait('L3', [0,0,1])

    full_vector_traits = fitch_with_trait_vector(Root, vector_traits) 

    assert(full_vector_traits.traits['root']==[0,0,0])
    assert(full_vector_traits.traits['I1']==[0,1,0])
