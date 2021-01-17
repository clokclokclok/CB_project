import time 
import progressbar as pb

bar = pb.ProgressBar(max_value=100)
for i in range(100):
    time.sleep(0.1)
    bar.update(i)

    