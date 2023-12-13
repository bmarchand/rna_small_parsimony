from rnadist.utils import C2_ILMedian

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
