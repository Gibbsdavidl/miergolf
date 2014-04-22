



    if s["printGraph"] == 1:
        x = coo_matrix(sparseLil)
        for i,j,v in zip(x.row, x.col, x.data):
            print "%s\t%s\t(%d, %d)\t%s" % (nodeDict[i], nodeDict[j], i,j,v)  
        exit(0)

