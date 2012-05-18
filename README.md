lsh
===

locality sensitive hashing

lsh is an indexing technique that makes it possible to search efficiently for nearest neighbours
amongst large collections of items, where each item is represented by a vector of some fixed dimension.
The algorithm is approximate but offers probabilistic guarantees i.e. with the right parameter
settings the results will rarely differ from doing a brute force search over your whole collection.
The search time will certainly be different though: lsh is useful because the complexity of lookups
becomes sublinear in the size of the collection.

In principle the algorithm is quite simple, but when I was getting to grips with it I couldn't find
any straightforward implementations just to see how it worked - so I wrote this one myself.  It's not
intended for use in production, but, depending on your requirements, you shouldn't find it too hard
to adapt it for production once you understand how it works.

The idea of lsh is to come up with a hashing scheme that maps closely neighbouring items to the same
bin, hence the "locality sensitive" part of its name.  The starting point is to pick a family of simple
hash functions.  Each member of this family is initialised with a different randomly chosen seed.  The
locality sensitive hash for an item is then constructed by joining together the values output by a
vector of k of these simple hash functions.  Intuitively you can see that once k is big enough, it's not
going to be a conincidence if two items map to the same sequence of k values.  On the other hand it's
possible that two items might be close neighbours but might still differ in one or two of the k values.
To increase recall we build L similar hash tables of this kind, where the items are mapped to bins using
k different hash functions for each table.  The lsh index is this set of L tables.  At search time we
apply all the simple hash functions to our query item and find its corresponding bin in each table.
Then we just have to check the exact distance of the items in those bins to decide which are the very
nearest neighbours.

The theory is explained quite clearly in this paper by the original creators of lsh:
A. Andoni and P. Indyk, "Near-optimal hashing algorithms for approximate nearest neighbor in high dimensions"
http://people.csail.mit.edu/indyk/p117-andoni.pdf

The main theoretical trick is to come up with a suitable family of simple hashes that gives the desired
result for the particular distance measure you're using to decide how similar items are to one another.
I've implemented the families most commonly used for L1, L2 and cosine distance.  There's also a simple
family that works with jaccard distance, which makes sense if you're items are really sets: in that case
the overall algorithm is often known as "minhash".  I've found minhashing to be extremely useful when
doing realtime lookups in collections of several million high-dimensional items.

If you want to move on to versions that are a bit closer to production code then I recommend you check
out these projects:
* https://github.com/andrewclegg/sketchy/blob/master/sketchy.py  (sketches i.e. lsh for cosine distance)
* https://github.com/hamilton/LocalitySensitiveHashing (minhash i.e. lsh for jaccard distance)
