def replace_line(filename, text):
    
    with open(filename) as file:
        lines = file.readlines()
    
    newlines = []

    for i in range(len(lines)):
        line = lines[i]
        if line[-2] == '0':
            lines[i] = lines[i].replace(" 0", "")
            newlines.append(lines[i])
    
    with open("2019Xlist_RxOff.txt", "w") as file:
        for line in newlines:
            file.write(line)


filename = "2019Xlist.txt"
text = " 1"

replace_line(filename, text)