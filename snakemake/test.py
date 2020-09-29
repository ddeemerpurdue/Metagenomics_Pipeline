mylist = ['dane', 'amanda', 'dillon']


for cnt, val in enumerate(mylist):
    if val == 'dane':
        newval = 'daners'
        mylist.pop(cnt)
        mylist.append(newval)

print(mylist)
