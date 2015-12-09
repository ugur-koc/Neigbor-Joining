#!/opt/local/bin/python2.7

from zss import simple_distance, Node

A = (
     Node("f")
     .addkid(Node("a")
             .addkid(Node("h"))
             .addkid(Node("c")
                     .addkid(Node("l"))))
     .addkid(Node("e"))
     )
B = (
     Node("f")
     .addkid(Node("a")
             .addkid(Node("d"))
             .addkid(Node("c")
                     .addkid(Node("b"))))
     .addkid(Node("e"))
     )


print simple_distance(A, B)