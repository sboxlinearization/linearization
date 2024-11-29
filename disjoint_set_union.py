class DisjointSetUnion:
    def __init__(self, init=None):
        self._parent = {}
        self._size = {}
        self._components = set()
        if init:
            for a, b in init:
                self.union(a, b)

    def find(self, a):
        while True:
            b = self._parent.get(a, a)
            if b == a:
                break
            self._parent[a] = self._parent.get(b, b)
            a = b
        return a

    def size(self, a):
        return self._size.get(self.find(a), 1)

    def union(self, a, b):
        a = self.find(a)
        b = self.find(b)

        if a == b:
            self._components.add(a)
            return

        sz_a = self.size(a)
        sz_b = self.size(b)
        if sz_a < sz_b:
            a, b = b, a

        self._components.discard(b)
        self._components.add(a)
        self._parent[b] = a
        self._size[a] = sz_a + sz_b


    def list(self):
        sets = {}
        for a in self._parent:
            root = self.find(a)
            if root not in sets:
                sets[root] = [root]
            if a != root:
                sets[root].append(a)
        return list(sets.values())

    def components(self):
        return list(self._components)


class RangeDSU:
    def __init__(self, n):
        self.n = n
        self._parent = [i for i in range(n)]
        self._size = {i:1 for i in range(n)}
        self._components = set(range(n))

    def find(self, a):
        while True:
            b = self._parent[a]
            if b == a:
                break
            self._parent[a] = self._parent[b]
            a = b
        return a

    def size(self, a):
        return self._size.get(self.find(a), 1)

    def union(self, a, b):
        a = self.find(a)
        b = self.find(b)

        if a == b:
            # already merged
            self._components.add(a)
            return

        sz_a = self._size[a]
        sz_b = self._size[b]
        if sz_a < sz_b:
            a, b = b, a

        self._components.discard(b)
        self._components.add(a)
        self._parent[b] = a
        del self._size[b]
        self._size[a] = sz_a + sz_b

    def list(self):
        sets = {}
        for a in range(self.n):
            root = self.find(a)
            if root not in sets:
                sets[root] = [a]
            else:
                sets[root].append(a)
        return list(sets.values())

    def components(self):
        return list(self._components)

class RangeDSU_with_min:
    def __init__(self, n):
        self.n = n
        self._parent = [i for i in range(n)]
        self._size = {i:1 for i in range(n)}
        self._components = set(range(n))
        self._min = {i:i for i in range(n)}
        self._max = {i:i for i in range(n)}

    def find(self, a):
        while True:
            b = self._parent[a]
            if b == a:
                break
            self._parent[a] = self._parent[b]
            a = b
        return a

    def size(self, a):
        return self._size.get(self.find(a), 1)

    def union(self, a, b):
        a = self.find(a)
        b = self.find(b)

        if a == b:
            self._components.add(a)
            return

        sz_a = self._size[a]
        sz_b = self._size[b]
        if sz_a < sz_b:
            a, b = b, a

        min_a = self._min[a]
        min_b = self._min[b]
        max_a = self._max[a]
        max_b = self._max[b]

        self._components.discard(b)
        self._components.add(a)
        self._parent[b] = a
        del self._size[b]        
        self._size[a] = sz_a + sz_b
        del self._min[b]
        self._min[a] = min(min_a, min_b)
        del self._max[b]
        self._max[a] = max(max_a, max_b) 

    def list(self):
        sets = {}
        for a in range(self.n):
            root = self.find(a)
            if root not in sets:
                sets[root] = [a]
            else:
                sets[root].append(a)
        return list(sets.values())

    def components(self):
        return list(self._components)

    def min(self,a):
        return self._min.get(self.find(a), a)

    def max(self,a):
        return self._max.get(self.find(a), a)

    def ec(self,a):
        root = self.find(a)
        cl = [i for i in range(self.n) if self.find(i) == root]
        return cl
        

