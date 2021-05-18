#! /usr/bin/python

import sys



def main():
    d = {}
    i = 1
    while 1:
        line = sys.stdin.readline()
        if line != '':
            if line[0] == 'v':
                xyz = list(map(float, line[2:].split(' ')))
                d[i] = [xyz[0], xyz[1], xyz[2]]
                # print(i, "--", xyz)
                i += 1
            
            elif line[0] == 'b': 
                pass

            elif line[0] == 'x': 
                val = list(map(int, line[2:].split(' ')))
                # print("del", val)
                del d[val[0]]

            elif line[0] == 'f': 
                ids = list(map(int, line[2:].split(' ')))
                # print(ids)
                for j in range(3):
                    if ids[j] not in d:
                        print ("oupsie", ids[j])
        else:
            if len(d) > 0:
                print("some vertices not finalised: ", d)
            else:
                print("all good 👍")
            break

if __name__ == "__main__":
    main()   
