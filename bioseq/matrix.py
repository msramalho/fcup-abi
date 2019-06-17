from .utils import module_exists


class Matrix:
    """Wrapper class to provide a numerical matrix interface with relevant and reusable functions"""

    def __init__(self, rows, cols=None, value=0, row_names=None, col_names=None):
        assert rows > 0, "The number of rows must be positive"
        if not cols: cols = rows  # square matrix
        assert cols > 0, "The number of cols must be positive"
        self.matrix = [[value for _ in range(cols)] for _ in range(rows)]
        assert not row_names or len(row_names) == rows, "The row names do no match the number of rows"
        assert not col_names or len(col_names) == cols, "The col names do no match the number of cols"
        self.rows = row_names or [""] * rows
        self.cols = col_names or [""] * cols

    def sum(self):
        """calculate the sum of the values in the matrix"""
        return sum(c for r in self for c in r)

    def max(self):
        """Maximum value in the matrix"""
        # return max(c for r in self for c in r)
        ma = max(c for r in self for c in r)
        return [(c, i, j) for i, r in enumerate(self) for j, c in enumerate(r) if c == ma]

    def min(self, ignore=None):
        """Minimum value in the matrix"""
        mi = min(c for r in self for c in r if c != ignore)
        return [(c, i, j) for i, r in enumerate(self) for j, c in enumerate(r) if c == mi]

    def square(self):
        """Test if matrix is squared"""
        return len(self) == len(self[0])

    def last(self):
        """Return last cell of matrix"""
        return self[-1][-1]

    def add_val(self, val):
        """add a value to all cells in the matrix. Returns self"""
        self.apply(lambda x, _i, _j: x + val)
        return self

    def mul_val(self, val):
        """multiply a value by all cells in the matrix. Returns self"""
        self.apply(lambda x, _i, _j: x * val)
        return self

    def symmetric(self):
        """Duplicate the values in the primary triangular matrix into the secondary one"""
        self.apply(lambda v, i, j: self[j][i] if j < i else v)
        return self

    def apply(self, operation):
        """general private method to apply an operation to each cell. The operation should receive (cell value, row index, col index)"""
        for i in range(len(self)):
            for j in range(len(self[i])):
                self[i][j] = operation(self[i][j], i, j)
        return self

    def display(self, save_to=False):  # pragma: no cover
        """Display the matrix in a plot if matplotlib is installed. If save_to is used the plot is saved and not shown."""
        if not module_exists("matplotlib"):
            return False
        import matplotlib.pyplot as plt
        plt.spy(self)
        if save_to:
            plt.savefig(save_to)
        else:
            plt.show()
        return True

    def graph(self, cut): # pragma: no cover
        """Draw a graph from the matrix, using only the superior triangle"""
        import networkx as nx
        g = nx.Graph()
        for i in range(len(self)): g.add_node(i)
        
        for i in range(len(self)):
            for j in range(i+1, len(self[0])):
                if self[i][j] < cut:
                    g.add_edge(i, j)
        
        nx.draw(g, with_labels=True)

    def __len__(self):
        """Get the length of the matrix - number of rows"""
        return len(self.matrix)

    def __setitem__(self, row, val):
        """set value of row"""
        assert len(val) == len(self.matrix[row]), "New row dimensions must match"
        self.matrix[row] = val

    def set_col(self, col, val):
        """set value of col"""
        assert len(val) == len(self.matrix), "New col dimensions must match"
        for i, row in enumerate(self): row[col] = val[i]

    def __getitem__(self, row):
        """get specific row of matrix"""
        return self.matrix[row]

    def __iter__(self):
        """iterator for the matrix, returns each row"""
        yield from self.matrix

    def __repr__(self):
        """for representation of Matrix"""
        return self.__str__()

    def __str__(self):
        """pretty print matrix"""
        # determine the necessary width
        use_r = bool(len(self.rows[0]))
        wl = max(len(str(c)) for r in self for c in r)
        wrl = min(max(map(len, self.rows)), 8)
        wcl = min(max(map(len, self.cols)), 8)
        wl = max(wl, wcl)
        wcl = wl + 1
        w = "%%%ds " % wl
        wr = "%%%ds:" % wrl if use_r else "%0s"
        wc = "%%%ds " % wcl
        res = " "*(wrl+1) if use_r else ""
        if use_r: res += w * len(self.cols) % tuple(map(lambda x: x[:wl], self.cols)) + "\n"
        for i, r in enumerate(self):
            res += wr % self.rows[i][:wl] + (w * len(r) % tuple(r)) + "\n"
        return res
