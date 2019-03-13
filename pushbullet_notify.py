from pushbullet import Pushbullet
import os
from sys import argv

pb = Pushbullet(os.environ["pushbullet_token"])
pb.push_link("FCUP-ABI pipeline %s" % argv[1], "https://github.com/msramalho/fcup-abi")