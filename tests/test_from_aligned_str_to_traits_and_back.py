from rnadist.utils import trait_vector_to_str, aligned_gapless_structures_to_trait_vectors

def test1():

    str_dict = {'1':'..((..))..',
                '2':'.(.(..))..',
                '3':'..((..).).'}

    trait_vectors, ordered_clade_list = aligned_gapless_structures_to_trait_vectors(str_dict)

    clade_list = [[3,4,5,6,7,8], 
                          [4,5,6,7], 
                          [2,3,4,5,6,7,8], 
                          [3,4,5,6,7,8,9]]

    clade_list = [frozenset(s) for s in clade_list]

    assert(set(ordered_clade_list)==set(clade_list))

    assert(trait_vector_to_str(trait_vectors.traits['1'], ordered_clade_list, 10)=='..((..))..')
    assert(trait_vector_to_str(trait_vectors.traits['2'], ordered_clade_list, 10)=='.(.(..))..')
    assert(trait_vector_to_str(trait_vectors.traits['3'], ordered_clade_list, 10)=='..((..).).')

def test2():

    str_dict = {'1' : '...................((((((....((((...((((...................))))(((((((((((...)))))))))))))))......))))))',
                '2' : '.....((((.....))))............((((((((((......))))))).)))......((((((((((.....))))))))))................',
                '3' : '.....((((.....)))).((((((....(((((((((((......)))))))..........(((((((((((...)))))))))))))))......))))))'}

    N = len(str_dict['1'])
    
    trait_vectors, ordered_clade_list = aligned_gapless_structures_to_trait_vectors(str_dict)
    
    assert(trait_vector_to_str(trait_vectors.traits['1'], ordered_clade_list, N)==str_dict['1'])
    assert(trait_vector_to_str(trait_vectors.traits['2'], ordered_clade_list, N)==str_dict['2'])
    assert(trait_vector_to_str(trait_vectors.traits['3'], ordered_clade_list, N)==str_dict['3'])
