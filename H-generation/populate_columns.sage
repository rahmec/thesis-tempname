from random import sample
import random

def generate_matrix(nrows, ncols, row_weight, column_weight):
    if nrows * row_weight != ncols * column_weight:
        print("Incompatible weights")
        return -1
    H = Matrix(nrows, ncols)
    available_indices = [0 for _ in range(ncols * column_weight)]

    for i in range(ncols * column_weight):
        available_indices[i] = i % nrows

    for col in range(ncols):
        rows = []
        while len(rows) < column_weight:
            row = random.choice(available_indices)
            while row in rows:
                available_indices2 = [value for value in available_indices if value not in rows]
                if len(available_indices2) == 0:
                    row = randint(0, nrows-1)
                    available_indices.remove(row)
                else:
                    row = random.choice(available_indices2)
                    available_indices.remove(row)
            rows.append(row)

        for row in rows:
            H[row, col] = 1

    print(H)

# THE MULTIPLICATION OF THE NUMBER OF ONES PER ROW MULTIPLIED BY THE NUMBER OF ROWSS
# MUST BE EQUAL TO
# THE MULTIPLICATION OF THE NUMBER OF ONES PER COLUMN MULTIPLIED BY THE NUMBER OF COLUMN
# => THERE IS A LIMITED NUMBER OF ALLOWED VALUES
# => PICKING A ROW VALUE IMPLIES A COLUMN VALUE
# => THEY DEPEND ON ONE ANOTHER
# Colweight = row_weight * nrows/ncols
# Se le colonne e le righe sono ugali stiamo apposto
# Se le colonne e le righe sono diverse dobbiamo scegliere un row_weight multiplo di ncols

generate_matrix(6,8,4,3)


