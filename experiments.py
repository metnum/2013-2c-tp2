# -*- coding: utf-8 -*-
import os
import collections
from subprocess import check_output
import numpy as np
from tp2 import f
from matplotlib import pyplot as plt


def force_position(pos, length):
    """
    Return the name of a force given its position
    """
    if pos == 0:
        return "H0"
    if pos == 1:
        return "V0"
    if pos < 4 * length - 1:
        return "F" + str(pos - 1)
    else:
        print pos
        return "V1"

def eq_position(row, n):
    if row % 2 == 1:
        vert_num = (row - 1) / 2
        return str(vert_num) + " H"
    else:
        vert_num = row /2
        return str(vert_num) + " V"


def check_matrix(m):
    # Horizontal and vertical link columns should add to 2
    # Diagonal links should add to 2/sqrt(2)

    def print_res(r):
        for result in r:
            print result[0], result[1]

    print "Getting sums of lower horizonal links"
    n = np.shape(m)[0] / 4
    sums = []
    for force in range(4, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:, f(force)]]))
    print_res(zip(range(4, 4*n, 4), sums))


    print "Getting sums of upper horizonal links"
    sums = []
    for force in range(6, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:, f(force)]]))
    print_res(zip(range(6, 4*n, 4), sums))


    print "Getting sums of vertical links"
    sums = []
    for force in range(3, 4*n -3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:, f(force)]]))
    print_res(zip(range(3, 4*n, 4), sums))

    print "Getting sums of diagonal links"
    sums = []
    for force in range(5, 4*n - 3, 4):
        sums.append(sum([abs(elem[0, 0]) for elem in m[:, f(force)]]))
    print_res(zip(range(5, 4*n, 4), sums))

    print "Showing links"
    for row in xrange(m.shape[0]):
        for column in xrange(m.shape[1]):
            if m[row, column] != 0.0:
                print ("Vertex " + str(row / 2) + " " + ("horizontal" if row % 2 == 0 else "vertical") + ": " +
                    force_position(column, m.shape[1]) + ", " + str(m[row, column]))


def format_result(result):
    """
    Given a result column vector, show the values assigned
    """
    if isinstance(result, list):
        n = len(result)
        for pos in range(0, n):
            print force_position(pos, n) + ": " + str(result[pos])
    else:
        n = result.shape[0]
        for pos in range(0, n):
            print force_position(pos, n) + ": " + str(result[pos, 0])


class BaseExperimento():
    # Executable args: span h n cargas...
    experiments = [] # List of experiments to run

    def max_force(self, prog_args, output):
        lines = [abs(float(line)) for line in output.split('\n')[1: -1]]
        return max(lines)

    def print_forces(self, prog_args, output):
        lines = [float(line) for line in output.split('\n')[1: -1]]
        span, h, n = prog_args[1: 4]
        C = prog_args[4:]
        format_result(lines)

    def __init__(self, experiments):
        self.experiments = experiments
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")
        self.resultados = []

        for args in self.experiments:
            prog_args = ["%s" % arg for arg in [executable] + args]
            output = check_output(prog_args)
            self.resultados.append(output)
            self.print_forces(prog_args, output)


class SpanStudyExp1(BaseExperimento):
    resultados = []
    C = [5, 5, 5, 5, 5, 5, 5]
    n = 8
    h = 3

    def max_force(self, prog_args, output):
        return float(output.split()[0])

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h
        spans = [float(v)/2 for v in xrange(self.n, 2 * 9 * self.n * self.h, self.n)]

        print ("Span, section width, section width / h, Max force")
        for span in spans:
            prog_args = ["%s" % arg for arg in [executable] + [
                span, self.h, self.n] + self.C]
            output = check_output(prog_args)
            self.resultados.append(output)
            print "%s, %s, %s, %s" % (span, (span/self.n), ((span/self.n)/self.h), self.max_force(prog_args, output))


class SpanStudyExp2(BaseExperimento):
    resultados = []
    C = [5] * 19
    n = 20
    h = 3

    def max_force(self, prog_args, output):
        return float(output.split()[0])

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h
        spans = [float(v)/2 for v in xrange(self.n, 2 * 9 * self.n * self.h, self.n)]

        print ("Span, section width, section width / h, Max force Link")
        for span in spans:
            prog_args = ["%s" % arg for arg in [executable] + [
                span, self.h, self.n] + self.C]
            output = check_output(prog_args)
            self.resultados.append(output)
            print "%s, %s, %s, %s" % (span, (span/self.n), ((span/self.n)/self.h), self.max_force(prog_args, output))


class AsymmetricWeightStudy(BaseExperimento):
    resultados = []
    C = ([0] * 3) + [47, 47] + [0] * 14
    n = 20
    h = 3

    def max_force(self, prog_args, output):
        return float(output.split()[0])

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h
        spans = [float(v)/2 for v in xrange(self.n, 2 * 9 * self.n * self.h, self.n)]

        print ("Span, section width, section width / h, Max force Link")
        for span in spans:
            prog_args = ["%s" % arg for arg in [executable] + [
                span, self.h, self.n] + self.C]
            output = check_output(prog_args)
            self.resultados.append(output)
            print "%s, %s, %s, %s" % (
                span, (span/self.n), ((span/self.n)/self.h), self.max_force(prog_args, output))



class SpanHistogramStudy(BaseExperimento):
    from matplotlib import pyplot as plt
    resultados = []
    C = [5] * 19
    n = 20
    h = 3

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h

        span = 200
        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', span, self.h, self.n] + self.C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s" % (self.max_force(prog_args, output))

        output = output.split('\n')
        cost = output[0]
        forces = output[1:]
        lines = [float(line) for line in forces if line]
        bins = range(-800, 800, 50)
        n, bins, patches = plt.hist(lines, bins, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=200, n=20, h=3, Ci = 5")
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_200.png')

        plt.show()
        span = 400

        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', span, self.h, self.n] + self.C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s" % (self.max_force(prog_args, output))

        output = output.split('\n')
        cost = output[0]
        forces = output[1:]
        lines = [float(line) for line in forces if line]
        bins = range(-1600,1600, 100)
        n, bins, patches = plt.hist(lines, bins, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=400, n=20, h=3, Ci = 5")
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_400.png')

        plt.show()


class SpanCentralWeightStudy(BaseExperimento):
    resultados = []
    C = ([0] * 9) + [5 * 19] + ([0] * 9)
    n = 20
    h = 3

    def max_force(self, prog_args, output):
        return float(output.split()[0])

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h
        spans = [float(v)/2 for v in xrange(self.n, 2 * 15 * self.n * self.h, self.n)]

        print ("Span, section width, section width / h, Max force Link")
        for span in spans:
            prog_args = ["%s" % arg for arg in [executable] + [
                span, self.h, self.n] + self.C]
            output = check_output(prog_args)
            self.resultados.append(output)
            print "%s, %s, %s, %s " % (
                span, (span/self.n), ((span/self.n)/self.h), self.max_force(prog_args, output))


class NHistStudy(BaseExperimento):
    resultados = []
    span = 100.0
    h = 3.0

    def max_force(self, prog_args, output):
        lines = [abs(float(line)) for line in output.split('\n')[1: -1] if line]
        return max(lines)

    def max_force_name(self, prog_args, output):
        lines = [abs(float(line)) for line in output.split('\n')[1: -1]]
        max_force = 0
        max_pos = -1
        for lineno, line in enumerate(lines):
            if line > max_force:
                max_force = line
                max_pos = lineno

        return force_position(max_pos, len(lines))

    def __init__(self):
        executable = os.path.join(os.getcwd(), "codigo/bin/tp2")

        # from 0.5 to 2 * 9 * n * h
        ns = [10, 100]

        print ("Span, section width, section width / h, Max force Link")
        n = 10
        total_c = 1000.0
        C = [total_c/(n - 1)] * int(n - 1)
        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', self.span, self.h, n] + C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s, %s, %s, %s, %s" % (
            self.span, (self.span/n), ((self.span/n)/self.h), self.max_force(prog_args, output),
            self.max_force_name(prog_args, output))

        output = output.split('\n')
        cost = output[0]
        forces = output[1:]
        lines = [float(line) for line in forces if line]

        bins = range(-1000, 1000, 200)
        x, bins, patches = plt.hist(lines, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=100, n=%s, h=3, Ci=%.2f" % (n, C[0]))
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_n' + unicode(n) + '_C1000.png')

        plt.show()

        total_c = 100.0

        C = [total_c/(n - 1)] * int(n - 1)
        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', self.span, self.h, n] + C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s, %s, %s, %s, %s" % (
            self.span, (self.span/n), ((self.span/n)/self.h), self.max_force(prog_args, output),
            self.max_force_name(prog_args, output))

        output = output.split('\n')
        cost = output[0]
        forces = output[1:]
        lines = [float(line) for line in forces if line]
        bins = range(-100, 100, 20)
        x, bins, patches = plt.hist(lines, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=100, n=%s, h=3, Ci=%.2f" % (n, C[0]))
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_n' + str(n) + '_C100.png')

        plt.show()

        n = 100

        total_c = 1000.0
        C = [total_c/(n - 1)] * int(n - 1)
        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', self.span, self.h, n] + C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s, %s, %s, %s, %s" % (
            self.span, (self.span/n), ((self.span/n)/self.h), self.max_force(prog_args, output),
            self.max_force_name(prog_args, output))

        lines = [float(line) for line in output.split('\n')[1: -1]]
        bins = range(-40000, 40000, 5000),
        x, bins, patches = plt.hist(lines, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=100, n=%s, h=3, Ci=%.2f" % (n, C[0]))
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_n' + unicode(n) + '_C1000.png')

        plt.show()

        total_c = 100.0

        C = [total_c/(n - 1)] * int(n - 1)
        prog_args = ["%s" % arg for arg in [executable] + [
            'display-forces', self.span, self.h, n] + C]
        output = check_output(prog_args)
        self.resultados.append(output)
        print "%s, %s, %s, %s, %s" % (
            self.span, (self.span/n), ((self.span/n)/self.h), self.max_force(prog_args, output),
            self.max_force_name(prog_args, output))

        lines = [float(line) for line in output.split('\n')[1: -1]]
        bins = range(-4000, 4000, 200)
        x, bins, patches = plt.hist(lines, histtype='bar')
        plt.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        plt.autoscale(True, 'both', False)
        plt.title(u"Distribución de fuerzas con Span=100, n=%s, h=3, Ci=%.2f" % (n, C[0]))
        plt.xlabel("Fuerza en viga")
        plt.ylabel("Cantidad")
        plt.grid(True)
        plt.legend()
        plt.savefig('informe/archivos/graficos/hist_n' + str(n) + '_C100.png')

        plt.show()

if __name__ == '__main__':
    # Fuerza maxima por span de puente
    # n=8, h=3, Ci=5
    #print "Study max span n=8, h=3, Ci=5"
    #study = SpanStudyExp1()

    #print "Study max span n=20, h=3, Ci=5"

    # Fueza máxima por span de puente
    # n=20, h=3, Ci=5
    #study = SpanStudyExp2()

    # Hist distribución de fuerzas span=200,n=20, h=3, Ci=5,
    study = SpanHistogramStudy()

    # Fuerza máxima carga central, span variable
    #study = SpanCentralWeightStudy()

    study = NHistStudy()
