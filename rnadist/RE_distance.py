#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

``RE_distance.py`` **module description**:


This module defines the function 'RE_distance(str1,str2)' 
to compute the RE distance between two RNA trees represented
as two strings str1 and str2 in dot-bracket notation.

The output of 'RE_distance(str1,str2)' is the RE distance.
"""


# arc substitution cost
def was(x,y,p,q):
    return abs(p-x)+abs(q-y)

# arc deletion cost
def wad(x,y):
    return 1

# arc insertion cost
def wai(p,q):
    return 1

# index for segment [x,y] of str1
# used for storage of distance in a 2D matrix
def ind1(n,x, y):
    if (y < x):
        return n*n
    else:
        return x*n+y

# index for segment [p,q] of str2
# used for storage of distance in a 2D matrix
def ind2(m,p, q):
    if (q < p):
        return m*m
    else:
        return p*m+q



def RE_distance(str1,str2):
    # lengths of strings
    n = len(str1)
    m=len(str2)

    #P1[i]==j means that (i,j) or (j,i) are paired in str1
    #P1[i]==-1 means that i is unpaired in str1
    P1=[]
    LP1 = [] #stack used to fill P1
    for i in range(0,n):
        P1.append(-1)
        LP1.append(-1)
        
    j = 0 #next position to add in stack
    k = 0 #top of stack
    for i in range(0,n):
        if(str1[i] == '('):
            LP1[j]=i
            k = j
            j+=1
        if(str1[i] == ')'):
            P1[LP1[k]]= i
            P1[i]= LP1[k]
            k-=1
            j-=1

    #P2[i]==j means that (i,j) or (j,i) are paired in str2
    #P2[i]==-1 means that i is unpaired in str2
    P2=[]
    LP2 = []#stack used to fill P2
    for i in range(0,m):
        P2.append(-1)
        LP2.append(-1)
        
    j = 0 #next position to add in stack
    k = 0 #top of stack
    for i in range(0,m):
        if(str2[i] == '('):
            LP2[j]=i
            k = j
            j+=1
        if(str2[i] == ')'):
            P2[LP2[k]]= i
            P2[i]= LP2[k]
            k-=1
            j-=1

    # I1[i][j]==1 means that (i,j) is an indexing pair of str1
    # indexing pair : a segment that do not cross any pair of str1
    # I1[i][j]== 0 means that (i,j) is not an indexing pair 
    I1=[]
    for i in range(0,n):
        I1.append([])
        for j in range(0,n):
            I1[i].append(0)
        j = i
        while(j < n):
            if(P1[j]==-1):
                I1[i][j]=1
                j += 1
            elif(P1[j] < i):
                j = n
            else:
                I1[i][P1[j]]=1
                j = P1[j] +1

    # I2[i][j]==1 means that (i,j) is an indexing pair of str2
    # indexing pair : a segment that do not cross any pair of str2
    # I2[i][j]== 0 means that (i,j) is not an indexing pair  
    I2=[]
    for i in range(0,m):
        I2.append([])
        for j in range(0,m):
            I2[i].append(0)
        j = i
        while(j < m):
            if(P2[j]==-1):
                I2[i][j]=1
                j += 1
            elif(P2[j] < i):
                j = m
            else:
                I2[i][P2[j]]=1
                j = P2[j] +1

    #Dynamic programming table indexed by pairs of indexing pairs
    D = []
    for i in range(n*n+1):
        D.append([])
        for j in range(m*m+1):
            D[i].append(0)

    #Initialisation
    D[n*n][m*m]=0

    for i in range(n):
        for x in range(n):
            y = x+i
            if(y < n):
                if(I1[x][y]==1):
                    if(P1[x] != -1):
                        D[ind1(n,x,y)][m*m]= wad(x,P1[x])+ D[ind1(n,x+1,P1[x]-1)][m*m] + D[ind1(n,P1[x]+1,y)][m*m]
                    else:
                        D[ind1(n,x,y)][m*m]=  D[ind1(n,x+1,y)][m*m]
                    
    
    for j in range(m):
        for p in range(m):
            q = p+j
            if(q < m):
                if(I2[p][q]==1):
                    if(P2[p] != -1):
                        D[n*n][ind2(m,p,q)]= wai(p,P2[p])+ D[n*n][ind2(m,p+1,P2[p]-1)] + D[n*n][ind2(m,P2[p]+1,q)]
                    else:
                        D[n*n][ind2(m,p,q)]=  D[n*n][ind2(m,p+1,q)]

    # Main recurrence formula
    for ii in range(n):
        for x in range(n):
            y = x+ii
            if(y < n):
                if(I1[x][y]==1):
                    for jj in range(m):
                        for p in range(m):
                            q = p+jj
                            if(q < m):
                                if(I2[p][q]==1):
                                    c = 0
                                    # if [x,y] starts with an unpaired base
                                    if(P1[x]==-1): 
                                        c = D[ind1(n,x+1,y)][ind2(m,p,q)]
                                    # if [p,q] starts with an unpaired base
                                    elif (P2[p]==-1):
                                        c = D[ind1(n,x,y)][ind2(m,p+1,q)]
                                    # if [x,y] and [p,q] starts with an arc
                                    else:
                                        # case where [x,y] and [p,q] are aligned
                                        c = was(x,P1[x],p,P2[p]) + D[ind1(n,x+1,P1[x]-1)][ind2(m,p+1,P2[p]-1)] + D[ind1(n,P1[x]+1,y)][ind2(m,P2[p]+1,q)]
                                        # case where [x,y] is deleted
                                        for i in range(p-1,q):
                                            if((i==p-1) or (i >= p and I2[p][i]==1)):
                                                tmp = wad(x,P1[x]) + D[ind1(n,x+1,P1[x]-1)][ind2(m,p,i)] + D[ind1(n,P1[x]+1,y)][ind2(m,i+1,q)]
                                                if(tmp < c):
                                                    c = tmp
                                        # case where [p,q] is deleted
                                        for i in range(x-1,y):
                                            if((i==x-1) or (i >= x and I1[x][i]==1)):
                                                tmp = wai(p,P2[p]) + D[ind1(n,x,i)][ind2(m,p+1,P2[p]-1)] + D[ind1(n,i+1,y)][ind2(m,P2[p]+1,q)]
                                                if(tmp < c):
                                                    c = tmp
                                    D[ind1(n,x,y)][ind2(m,p,q)] =c
                                    
    return D[ind1(n,0,n-1)][ind2(m,0,m-1)]

  

