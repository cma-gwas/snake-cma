# Connect to task counter server
import zmq
import argparse

def counter_next(job_id):
    socket.send(bytes(job_id, encoding="utf-8"))
    b_value = socket.recv()
    return int.from_bytes(b_value, "big")

def send_quit():
    socket.send(bytes("quit", encoding="utf-8"))
    

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--port", type=int, default=5555, help="Server port")
args = parser.parse_args()

context = zmq.Context()
socket = context.socket(zmq.REQ)

print("Connecting to counter server on port {}...".format(args.port))
socket.connect("tcp://localhost:{}".format(args.port))

#resp = counter_next("quit")
#print(resp)
send_quit()
