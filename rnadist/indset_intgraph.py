import operator
import sys




#returns the index i of the sorted_list such that sorted_list[i] < x
#Takes linear time, pretty bad.  Could be optimized to log n with binary search
def find_largest_index_under( sorted_list, x ):
    
    cur = None
    
    for i in range(len(sorted_list)):
        if sorted_list[i] < x:
            cur = i
        else:
            return cur

    return cur
   

def build_solution( B ):
    
    sol = list()
    cur = len(B) - 1 
    
    while cur != None:
        if B[cur][0] != None:
            sol.append(B[cur][0])
        cur = B[cur][1]
        
    return sol
        
   
    



# intervals : list of triples of the form (a, b, w) where a = left endpoint, b = right endpoint, w = weight
# return_weight_only : if True, returns only the weight of the max indset, if False, returns an actual solution
# already_sorted : set to True if intervals is sorted by right endpoints (in which case the sorting will be skipped), and False otherwise
def get_max_interval_indset( intervals, return_weight_only = False, already_sorted = False ):

    #sort according to right endpoint
    intervals.sort(key=operator.itemgetter(1))

    #ends = sorted list of all possible right endpoints - we will iterate over these
    ends = []
    for i in range(len(intervals)):
        if len(ends) == 0 or ends[-1] != intervals[i][1]:
            ends.append(intervals[i][1])

    maxend = ends[-1]

    #debug stuff
    #print(intervals)
    #print(ends)
    


    #for each possible e_index in [0, len(ends) - 1], V[e_index] is the size of max indset of intervals entirely contained in [1..ends[e_index] ] 
    V = [0] * (len(ends))
    
    #B is the backtracking array, B[e_index] = ( I, e2 ) with I is the interval to include whose end is ends[e_index] (or None) 
    #and e2 the previous end index to go to
    B = [None] * (len(ends))

    interval_index = 0

    for (e_index, e) in enumerate(ends):
        
        #print( "ends[" + str(e_index) + "]=" + str(e))
        
        
        # case 1, best weight does not use anything that ends at e
        if e_index > 0:
            V[e_index] = V[ e_index - 1 ]
            B[ e_index ] = ( None, e_index - 1 )

            
        # case 2, best weight uses something that ends at e, we consider all intervals that end at e
        #         since intervals are sorted, those intervals are grouped together, the first being at interval_index (supposed to be, which we check)
        if intervals[ interval_index ][1] != e:
            print("Error, item " + str(interval_index) + " should have end " + str(e) + " but has " + str(intervals[ interval_index ][1]))
            sys.exit()
        
        
        while interval_index < len(intervals) and intervals[ interval_index ][1] == e:
            
            best_weight = intervals[ interval_index ][2]
            
            start = intervals[ interval_index ][0]
            prev_end_index = find_largest_index_under( ends, start )
            
            if prev_end_index != None:
                best_weight += V[prev_end_index]
                
            if best_weight > V[e_index]: 
                V[e_index] = best_weight
                B[ e_index ] = (intervals[ interval_index ], prev_end_index)
            
            interval_index += 1

    if return_weight_only:
        return V[-1]
    else:
        sol = build_solution( B ) 
        return sol
            

# format (start, end, weight)
intervals = [ (1, 3, 4), (5, 7, 5), (8, 9, 3), (6, 8, 2), (2, 7, 10) ]
sol = get_max_interval_indset( intervals )

weight = sum( i[2] for i in sol )

print( sol )
print(f"weight = {weight}")
