import subprocess

import numpy as np

r = 0x1
g = 0x1000
b = 0x100000
a = 0x10000000

class VideoSink(object) :

    def __init__( self, size, filename="output", rate=20, byteorder="rgba" ) :
            self.size = size
            cmdstring  = ('/usr/local/bin/mencoder',
                    '/dev/stdin',
                    '-demuxer', 'rawvideo',
                    '-rawvideo', 'w=%i:h=%i'%size[::-1]+":fps=%i:format=%s"%(rate,byteorder),
                    '-o', filename+'.avi',
                    '-ovc', 'lavc',
                    "-force-avi-aspect",'1.0',
                    '-sws','2',
                    '-vf','scale=800:800'
                    )
            self.p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE, shell=False)

    def run(self, image) :
            #assert image.shape == self.size
            self.p.stdin.write(image.tostring())
    def close(self) :
            self.p.stdin.close()


def numpy2rgb(data):

    out = np.uint32(r*data+g*data+b*data+255*a)
    return out

def rgba2hex(s):
    out = np.uint32(r*s[0]+g*s[1]+b*s[2]+a*s[3])
    return out
