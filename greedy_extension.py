from heapq import heappush, heappop
from random import choice, random as random_ratio, shuffle
from disjoint_set_union import RangeDSU

def GreedyRunHeap(Alg, init=(), force=False, debug=0, profile=False):
    n = Alg.n
    m = Alg.m
    S = Alg.S
    
    LS = {0}
    LSx = {0}
    
    if profile:
        profile = [1]

    cliques = RangeDSU(2**n)
    queue = [(-1, random_ratio(), i) for i in range(2**n)]
    shuffle(queue)
    
    ddt2 = Alg.compute_kddt(k=2)
    
    init = list(init)    
    
    for i in range(n):
        if debug: print("iteration", i, "\n========")

        cl1 = cl2 = None
        while len(init) >= 2:
            cl1 = cliques.find(init.pop())
            cl2 = cliques.find(init.pop())
            init.append(cl1)
            if cl1 == cl2:
                cl1 = cl2 = None
                continue
            break
        else:
            if init and force:
                cl1 = init[0]
            
        if cl1 is None:
            while True:
                sz1, _, cl1 = heappop(queue)
                sz1 = -sz1
                if cliques.size(cl1) == sz1:
                    break
        
        # find max-size clique from another coset
        if cl2 is None:
            tmp = []
            while True:
                sz2, _, cl2 = heappop(queue)
                sz2 = -sz2
                if cliques.size(cl2) != sz2:
                    continue
                tmp.append((-sz2, random_ratio(), cl2))
                dx = cl1 ^ cl2
                if dx not in LSx:
                    break         
            for q in tmp:
                heappush(queue, q)
        
        if profile:
            profile.append(sz1 + sz2)
    
        # cl identifier is one of the elements itself
        x1 = cl1
        x2 = cl2
        dxdy0 = ((x1 ^ x2) << m) | (S[x1] ^ S[x2])
        dx0 = x1 ^ x2

        # merge cliques
        for dxdy in LS:
            dxdy ^= dxdy0        
            for x1, x2 in ddt2[dxdy >> m][dxdy & (2**m-1)]:
                cliques.union(x1, x2)
                sz = cliques.size(x1)
                heappush(queue, (-sz, random_ratio(), x1))
              
        LS = LS | {v ^ dxdy0 for v in LS}
        LSx = LSx | {v ^ dx0 for v in LSx}
        
    cls = cliques.list()
    mx = max(len(cl) for cl in cls)
    for cl in cls:
        if len(cl) == mx:
            break
    else:
        assert 0
    if profile:
        return cl, profile
    return cl
