import sys 

allsh = "exec_all.sh"

def filename_dfs(filenames,n,i,namevars,selectedvars,contents,m,j,indice,p,k,contvars,name):
    if i == n:
        with open(name,mode="w") as f:
            content_dfs(selectedvars,contents,m,0,indice,p,0,contvars,contents[0],f)
        with open(allsh,mode="a") as f:
            f.write("./" + name+"\n")
    else:
        for var in namevars[i]:
            filename_dfs(filenames,n,i+1,namevars,selectedvars+" "+var,contents,m,j,indice,p,k,contvars,name+var+filenames[i+1])


def content_dfs(selectedvars,contents,m,j,indice,p,k,contvars,command,f):
    if j == m:
        f.write(command+"\n")
    else:
        # from selectedvars
        if indice[j] != -1:
            var = list(selectedvars.split())[indice[j]]
            content_dfs(selectedvars,contents,m,j+1,indice,p,k,contvars,command+var+contents[j+1],f)
        # original var
        else:
            for var in contvars[k]:
                content_dfs(selectedvars,contents,m,j+1,indice,p,k+1,contvars,command+var+contents[j+1],f)

if __name__ == "__main__":
    print("--- Shell script generator for computer experiments ---")
    print("generates exec_all.sh which executes all generated scripts")
    filename = input("base file name (type * for variables): ")

    filenames = list(filename.split("*"))
    n = len(filenames)-1

    namevars = []
    for i in range(n):
        var = input("{}-th variable list for filename seperated by space: ".format(i))
        namevars.append(list(var.split()))

    content = input("base command in the file content (type * for variables): ")
    contents = list(content.split("*"))
    m = len(contents)-1

    indice = list(map(int,input("index list for each * (i if it's same as i-th var in filename otherwise -1): ").split()))

    p = indice.count(-1)
    contvars = []
    for i in range(p):
        var = input("{}-th variable list for command seperated by space: ".format(i))
        contvars.append(list(var.split()))


    print("generating shell scripts...")
    filename_dfs(filenames,n,0,namevars,"",contents,m,0,indice,p,0,contvars,filenames[0])
    print("done.")





    
