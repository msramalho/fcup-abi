from .utils import module_exists


class Matrix:
    """Wrapper class to provide a numerical matrix interface with relevant and reusable functions"""

    def __init__(self, rows, cols, value=0):
        assert rows > 0, "The number of rows must be positive"
        assert cols > 0, "The number of cols must be positive"
        self.matrix = [[value for _ in range(cols)] for _ in range(rows)]

    def sum(self):
        """calculate the sum of the values in the matrix"""
        return sum(c for r in self for c in r)

    def max(self):
        """Maximum value in the matrix"""
        return max(c for r in self for c in r)

    def min(self):
        """Minimum value in the matrix"""
        return min(c for r in self for c in r)

    def square(self):
        """Test if matrix is squared"""
        return len(self) == len(self[0])

    def add_val(self, val):
        """add a value to all cells in the matrix. Returns self"""
        self.apply(lambda x, _i, _j: x + val)
        return self

    def mul_val(self, val):
        """multiply a value by all cells in the matrix. Returns self"""
        self.apply(lambda x, _i, _j: x * val)
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

    def __str__(self):
        """pretty print matrix"""
        # determine the necessary width
        w = max(len(str(c)) for r in self for c in r)
        res = ""
        for r in self: res += (("%%%ds " % w) * len(r) % tuple(r)) + "\n"
        return res
