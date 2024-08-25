import zmq
import atomics
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--port", type=int, default=5555, help="Listen port")
args = parser.parse_args()

context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:{}".format(args.port))

global_counter = atomics.atomic(width=4, atype=atomics.INT)

print("Counter server is running on port {} ...".format(args.port))

stopped = False

while not stopped:
    #  Wait for next request from client
    message = socket.recv()
    #print("Received request: %s" % message)
    if str(message, 'utf-8').lower() == "quit":
        stopped = True
        #socket.send((0).to_bytes(4, "big"))
    else:
        #  Do some 'work'
        next_value = global_counter.fetch_inc()
        print("Received request: %s, counter -> %s" % (message, next_value))

        #  Send reply back to client
        #socket.send(next_value)
        socket.send(next_value.to_bytes(4, "big"))

print("Counter server is stopped (port={})".format(args.port))
#context.destroy()
