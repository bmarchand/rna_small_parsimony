from rnadist.utils import C2_ILMedian, check_well_parenthesized

def test1_C2_ILMedian():

    s = '...((...))...'

    m = C2_ILMedian([s,s,s])

    assert(s==m)

def test2_C2_ILMedian():

    s1 = '...((...))...'
    s2 = '...(.(..))...'
    s3 = '...((...))...'

    m = C2_ILMedian([s1,s2,s3])

    assert(s1==m)

def test3_buggy_case():

    strs = ['((.)(.)(.))', 
            '((.)(.)(.))', 
            '((.)(.)(.))']
    
    m = C2_ILMedian(strs)

def test4():
    strs = ['.(((((......)))))(((((())))))', '.(((((......)))))(((((())))))', '.(((((......)))))(((((())))))']
    
    m = C2_ILMedian(strs)

    check_well_parenthesized(m)


def test5():
#                     1111111111
#           01234567890123456789
    strs = ['((((..(())))))....',
            '((((..(())....))))', 
            '.(((((.(.).)).))).']

    for s in strs:
        check_well_parenthesized(s)
    
    m = C2_ILMedian(strs)

    check_well_parenthesized(m)

def test6():
    #       0123456789
    strs = ['...(.)()', 
            '(.)(.)()', 
            '(.)(.)()']
    
    m = C2_ILMedian(strs)

    check_well_parenthesized(m)
