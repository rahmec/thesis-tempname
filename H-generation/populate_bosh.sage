from random import sample
import random

def generate_matrix(nrows, ncols, row_weight, column_weight):
    if nrows * row_weight != ncols * column_weight:
        print("Incompatible weights")
        return -1
    H = Matrix(nrows, ncols)
    rows_to_fill = {}
    cols_to_fill = {}
    pickable_cols = [ [x for x in range(ncols)] for _ in range(nrows)]

    for i in range(nrows):
        rows_to_fill[i] = row_weight

    for i in range(ncols):
        cols_to_fill[i] = column_weight

    while rows_to_fill and cols_to_fill:
        row = random.choice(list(rows_to_fill.keys()))
        col = random.choice([value for value in list(cols_to_fill.keys()) if value in pickable_cols[row]])
        pickable_cols[row].remove(col)

        H[row, col] = 1

        if rows_to_fill[row] != 1:
            rows_to_fill[row] = rows_to_fill[row]-1
        else:
            rows_to_fill.pop(row)

        if cols_to_fill[col] != 1:
            cols_to_fill[col] = cols_to_fill[col]-1
        else:
            cols_to_fill.pop(col)
            #print(pickable_cols)
            #print('removing col ', col)
            for i in range(len(pickable_cols)):
                if col in pickable_cols[i]:
                    pickable_cols[i].remove(col)
            #print(pickable_cols)
            #print(rows_to_fill)
            #print(cols_to_fill)
            #print(H)
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

generate_matrix(600,800,40,30)

