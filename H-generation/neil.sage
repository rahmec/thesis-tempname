import random

def generate_matrix(m, n, wr, wc):
    H = Matrix(m,n)
    total_ones_indices = [[] for _ in range(m)]
    for i in range(n):
        available_indices = [x for x in range(m)]
        for _ in range(wc):
            row = random.choice(available_indices)
            H[row, i]=1
            available_indices.remove(row)
            total_ones_indices[row].append(i)
    
    underweighted_rows = []
    overweighted_rows = []
    for i in range(m):
        #if len(total_ones_indices[i]) < 2:
        #    underweighted_rows.append(i)
        #    zero_indices = [r for r in range(n) if r not in total_ones_indices[i]]
        #    for _ in range(2-len(total_ones_indices[i])):
        #        col = random.choice(zero_indices)
        #        zero_indices.remove(col)
        #        H[i,col]=1
        #        total_ones_indices[i].append(col)
        if len(total_ones_indices[i]) < wr:
            underweighted_rows.append(i)
        elif len(total_ones_indices[i]) > wr:
            overweighted_rows.append(i)

    print(overweighted_rows)
    print(underweighted_rows)

    for ov_row in overweighted_rows:
        for _ in range(len(total_ones_indices[ov_row])-wr):
            found_flag = False
            available_indices = total_ones_indices[ov_row]
            while not found_flag:
                col = random.choice(available_indices)
                #print(H)
                for un_row in underweighted_rows:
                    if col not in total_ones_indices[un_row]:
                        H[un_row, col] = 1
                        H[ov_row, col] = 0
                        total_ones_indices[un_row].append(col)
                        total_ones_indices[ov_row].remove(col)
                        if len(total_ones_indices[un_row]) == wr:
                            underweighted_rows.remove(un_row)
                        found_flag = True
                        break
                if col in available_indices:
                    available_indices.remove(col)
    
    print(H)
    return H
                    

    
H = generate_matrix(800,1200,600,400) 
print(len([a for a in H.list() if a]))
for row in H:
    if len([a for a in row.list() if a]) != 600:
        print("UN PESO SBAGLIATO TROVATO!")
