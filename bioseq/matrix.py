
class Matrix:
    """Wrapper class to provide a numerical matrix interface with relevant and reusable functions"""

    def __init__(self, rows, cols, value=0):
        assert rows > 0, "The number of rows must be positive"
        assert cols > 0, "The number of cols must be positive"
        self.matrix = [[value for _ in range(cols)] for _ in range(rows)]

    def sum(self):
        """calculate the sum of the values in the matrix"""
        return sum(c for r in self for c in r)

    def add_val(self, val):
        """add a value to all cells in the matrix. Returns self"""
        self._apply(lambda x: x + val)
        return self

    def mul_val(self, val):
        """multiply a value by all cells in the matrix. Returns self"""
        self._apply(lambda x: x * val)
        return self

    def _apply(self, operation):
        """general private method to apply an operation to each cell"""
        for i in range(len(self)):
            for j in range(len(self[i])):
                self[i][j] = operation(self[i][j])
        return self

    def max(self):
        """Maximum value in the matrix"""
        return max(c for r in self for c in r)

    def min(self):
        """Minimum value in the matrix"""
        return min(c for r in self for c in r)

    def square(self):
        """Test if matrix is squared"""
        return len(self) == len(self[0])

    def __len__(self):
        """Get the length of the matrix - number of rows"""
        return len(self.matrix)

    def __getitem__(self, index):
        """get specific row of matrix"""
        return self.matrix[index]

    def __iter__(self):
        """iterator for the matrix, returns each row"""
        yield from self.matrix
