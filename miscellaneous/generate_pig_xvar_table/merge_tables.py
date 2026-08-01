X_nint = 5

filename = 'pig_table.dat'

with open(filename, "wb") as outfile:
    for i in range(X_nint + 1):
        fl = f"id_{i}/pig_table.dat"
        with open(fl, "rb") as infile:
            outfile.write(infile.read())
