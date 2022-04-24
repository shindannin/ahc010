#!/usr/bin/python3

import subprocess


file = open('scoresheet/scoresheet_now.txt', 'w')  #書き込みモードでオープン

# for testSeed in range(startSeed,startSeed+numSeeds):
#    res = subprocess.run(["java","-jar","./Tester.jar","-exec","visualstudio/Release/marathon_main.exe","-seed",str(testSeed),"-novis"], stdout=subprocess.PIPE)
#    text = res.stdout.decode('utf-8')
#    sp = text.split("Score = ")
#    score = sp[-1].strip()
#    totalScore += float(score)
#    print(testSeed, score)
#    file.write(str(sp[-1]))


import time
from multiprocessing import Pool

def nijou(case_id):
    # 実行してoutputを出力する場合

    print(res)
    return(res)

if __name__ == "__main__":

    numCases = 40
    pool = Pool(8)
    result = pool.map(nijou, range(numCases))

    for case_id in range(numCases):
        # inputとoutputの両方のファイルを
        res = subprocess.run(["python3","./tester/judge.py",f"./data/in/input_{case_id}.txt",f"./data/out/output{case_id}.txt"], stdout=subprocess.PIPE)


    print(result)
