#!/usr/bin/env python3

import argparse
from itertools import combinations_with_replacement, combinations


# Tree node class
class Node:
    # tpl_or_lab - tuple=internal node
    #           or string=a leaf (label)
    # parent - defines a parent of a non-root node
    def __init__(self, tpl_or_lab, parent=None):
        if type(tpl_or_lab) == tuple:
            self.left = Node(tpl_or_lab[0], self)
            self.right = Node(tpl_or_lab[1], self)
            self.label = None
        else:
            self.label = tpl_or_lab
        self.parent = parent
        self.src = tpl_or_lab

    # is a leaf
    def leaf(self):
        return self.label

    def __str__(self):
        if self.leaf():
            return str(self.label)
        return "(%s,%s)" % (str(self.left), str(self.right))


def shapes(n):
    if n == 1:
        return [Node('X')]
    elif n == 2:
        return [Node(('X', 'X'))]
    else:
        nodes = []
        if n % 2 > 0:
            for k in range(1, (n - 1) // 2 + 1):
                a = shapes(k)
                b = shapes(n - k)
                nodes_comb = combinations([*a, *b], 2)
                for nc in nodes_comb:
                    nod = Node(nc)
                    if str(nod).count(',') <= (n - 1):
                        nodes.append(nod)
        else:
            nodes = []
            for k in range(1, n // 2):
                a = shapes(k)
                b = shapes(n - k)
                nodes_comb = combinations([*a, *b], 2)
                for nc in nodes_comb:
                    nod = Node(nc)
                    if str(nod).count(',') <= (n - 1):
                        nodes.append(nod)

            a = shapes(n // 2)
            nodes_comb = combinations_with_replacement(a, 2)

            for nc in nodes_comb:
                nod = Node(nc)
                if str(nod).count(',') <= (n - 1):
                    nodes.append(nod)
        return nodes


def gen_shapes(n):
    print(f'Shapes for n = {n}\n')

    shapes_generated = shapes(n)
    results = set()
    for s in shapes_generated:
        if str(s).count(',') == (n - 1):
            shap = str(s).replace('X', '')
            if shap not in results:
                results.add(shap)
                print(shap)

    print(f'\nSet size = {len(results)}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tree Shape Generator')
    parser.add_argument('num', type=int, help='shapes for specific n')

    args = parser.parse_args()
    gen_shapes(args.num)
