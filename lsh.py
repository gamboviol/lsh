"""
LSH Locality Sensitive Hashing
- indexing for nearest neighbour searches in sublinear time

simple tutorial implementation based on
A. Andoni and P. Indyk, "Near-optimal hashing algorithms for approximate nearest neighbor in high dimensions"
http://people.csail.mit.edu/indyk/p117-andoni.pdf
"""

import random
from collections import defaultdict
from operator import itemgetter

class LSHIndex:

    def __init__(self,hash_family,k,L):
        self.hash_family = hash_family
        self.k = k
        self.L = 0
        self.hash_tables = []
        self.resize(L)

    def resize(self,L):
        """ update the number of hash tables to be used """
        if L < self.L:
            self.hash_tables = self.hash_tables[:L]
        else:
            # initialise a new hash table for each hash function
            hash_funcs = [[self.hash_family.create_hash_func() for h in xrange(self.k)] for l in xrange(self.L,L)]
            self.hash_tables.extend([(g,defaultdict(lambda:[])) for g in hash_funcs])

    def hash(self,g,p):
        return self.hash_family.combine([h.hash(p) for h in g])

    def index(self,points):
        """ index the supplied points """
        self.points = points
        for g,table in self.hash_tables:
            for ix,p in enumerate(self.points):
                table[self.hash(g,p)].append(ix)
        # reset stats
        self.tot_touched = 0
        self.num_queries = 0

    def query(self,q,metric,max_results):
        """ find the max_results closest indexed points to q according to the supplied metric """
        candidates = set()
        for g,table in self.hash_tables:
            matches = table.get(self.hash(g,q),[])
            candidates.update(matches)
        # update stats
        self.tot_touched += len(candidates)
        self.num_queries += 1
        # rerank candidates
        candidates = [(ix,metric(q,self.points[ix])) for ix in candidates]
        candidates.sort(key=itemgetter(1))
        return candidates[:max_results]

    def get_avg_touched(self):
        """ mean number of candidates inspected per query """
        return self.tot_touched/self.num_queries

class L1HashFamily:

    def __init__(self,w,d):
        self.w = w
        self.d = d

    def create_hash_func(self):
        # each L1Hash is initialised with a different random partition vector
        return L1Hash(self.rand_partition(),self.w)

    def rand_partition(self):
        return [random.uniform(0,self.w) for i in xrange(self.d)]

    def combine(self,hashes):
        """
        combine hash values naively with str()
        - in a real implementation we can mix the values so they map to integer keys
        into a conventional map table
        """
        return str(hashes)

class L1Hash:

    def __init__(self,S,w):
        self.S = S
        self.w = w

    def hash(self,vec):
        # use str() as a naive way of forming a single value
        return str([int((vec[i]-s)/self.w) for i,s in enumerate(self.S)])

def L1_norm(u,v):
        return sum(abs(u[i]-v[i]) for i in xrange(len(u)))

def dot(u,v):
    return sum(ux*vx for ux,vx in zip(u,v))

class L2HashFamily:

    def __init__(self,w,d):
        self.w = w
        self.d = d

    def create_hash_func(self):
        # each L2Hash is initialised with a different random projection vector and offset
        return L2Hash(self.rand_vec(),self.rand_offset(),self.w)

    def rand_vec(self):
        return [random.gauss(0,1) for i in xrange(self.d)]

    def rand_offset(self):
        return random.uniform(0,self.w)

    def combine(self,hashes):
        """
        combine hash values naively with str()
        - in a real implementation we can mix the values so they map to integer keys
        into a conventional map table
        """
        return str(hashes)

class L2Hash:

    def __init__(self,r,b,w):
        self.r = r
        self.b = b
        self.w = w

    def hash(self,vec):
        return int((dot(vec,self.r)+self.b)/self.w)

def L2_norm(u,v):
        return sum((ux-vx)**2 for ux,vx in zip(u,v))**0.5

class CosineHashFamily:

    def __init__(self,d):
        self.d = d

    def create_hash_func(self):
        # each CosineHash is initialised with a random projection vector
        return CosineHash(self.rand_vec())

    def rand_vec(self):
        return [random.gauss(0,1) for i in xrange(self.d)]

    def combine(self,hashes):
        """ combine by treating as a bitvector """
        return sum(2**i if h > 0 else 0 for i,h in enumerate(hashes))

class CosineHash:

    def __init__(self,r):
        self.r = r

    def hash(self,vec):
        return self.sgn(dot(vec,self.r))

    def sgn(self,x):
        return int(x>0)

def cosine_distance(u,v):
    return 1 - dot(u,v)/(dot(u,u)*dot(v,v))**0.5

class LSHTester:
    """
    grid search over LSH parameters, evaluating by finding the specified
    number of nearest neighbours for the supplied queries from the supplied points
    """

    def __init__(self,points,queries,num_neighbours):
        self.points = points
        self.queries = queries
        self.num_neighbours = num_neighbours

    def run(self,name,metric,hash_family,k_vals,L_vals):
        """
        name: name of test
        metric: distance metric for nearest neighbour computation
        hash_family: hash family for LSH
        k_vals: numbers of hashes to concatenate in each hash function to try in grid search
        L_vals: numbers of hash functions/tables to try in grid search
        """
        exact_hits = [[ix for ix,dist in self.linear(q,metric,self.num_neighbours+1)] for q in self.queries]

        print name
        print 'L\tk\tacc\ttouch'
        for k in k_vals:        # concatenating more hash functions increases selectivity
            lsh = LSHIndex(hash_family,k,0)
            for L in L_vals:    # using more hash tables increases recall
                lsh.resize(L)
                lsh.index(self.points)

                correct = 0
                for q,hits in zip(self.queries,exact_hits):
                    lsh_hits = [ix for ix,dist in lsh.query(q,metric,self.num_neighbours+1)]
                    if lsh_hits == hits:
                        correct += 1
                print "{0}\t{1}\t{2}\t{3}".format(L,k,float(correct)/100,float(lsh.get_avg_touched())/len(self.points))

    def linear(self,q,metric,max_results):
        """ brute force search by linear scan """
        candidates = [(ix,metric(q,p)) for ix,p in enumerate(self.points)]
        return sorted(candidates,key=itemgetter(1))[:max_results]


if __name__ == "__main__":

    # create a test dataset of vectors of non-negative integers
    d = 5
    xmax = 20
    num_points = 1000
    points = [[random.randint(0,xmax) for i in xrange(d)] for j in xrange(num_points)]

    # seed the dataset with a fixed number of nearest neighbours
    # within a given small "radius"
    num_neighbours = 2
    radius = 0.1
    for point in points[:num_points]:
        for i in xrange(num_neighbours):
            points.append([x+random.uniform(-radius,radius) for x in point])

    # test lsh versus brute force comparison by running a grid
    # search for the best lsh parameter values for each family
    tester = LSHTester(points,points[:num_points/10],num_neighbours)

    args = {'name':'L2',
            'metric':L2_norm,
            'hash_family':L2HashFamily(10*radius,d),
            'k_vals':[2,4,8],
            'L_vals':[2,4,8,16]}
    tester.run(**args)

    args = {'name':'L1',
            'metric':L1_norm,
            'hash_family':L1HashFamily(10*radius,d),
            'k_vals':[2,4,8],
            'L_vals':[2,4,8,16]}
    tester.run(**args)

    args = {'name':'cosine',
            'metric':cosine_distance,
            'hash_family':CosineHashFamily(d),
            'k_vals':[16,32,64],
            'L_vals':[2,4,8,16]}
    tester.run(**args)