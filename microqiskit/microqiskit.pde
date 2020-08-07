// This is the Processing version of Qiskit. For the full version, see qiskit.org.
// It has many more features, and access to real quantum computers.
import java.util.List;
import java.util.Arrays;
import java.util.Map;

class QuantumCircuit {
  int numQubits;
  int numClbits;
  List<List> data = new ArrayList();

  QuantumCircuit(int n, int m) {
    this.numQubits = n;
    this.numClbits = m;
  }

  QuantumCircuit(int n) {
    this(n, 0);
  }


  // Applies an x gate to the given qubit.
  QuantumCircuit x(int q) {
    this.data.add(Arrays.asList("x", new Integer(q)));
    return this;
  }


  // Applies an rx gate to the given qubit by the given angle.
  QuantumCircuit rx(double theta, int q) {
    this.data.add(Arrays.asList("rx", new Double(theta), new Integer(q)));
    return this;
  }

  // Applies an h gate to the given qubit.
  QuantumCircuit h(int q) {
    this.data.add(Arrays.asList("h", new Integer(q)));
    return this;
  }


  // Applies a cx gate to the given source and target qubits.
  QuantumCircuit cx(int s, int t) {
    this.data.add(Arrays.asList("cx", new Integer(s), new Integer(t)));
    return this;
  }


  // Applies an rz gate to the given qubit by the given angle.
  QuantumCircuit rz(double theta, int q) {
    // This gate is constructed from `h` and `rx`.
    this.h(q);
    this.rx(theta, q);
    this.h(q);
    return this;
  }


  // Applies an ry gate to the given qubit by the given angle.
  QuantumCircuit ry(double theta, int q) {
    // This gate is constructed from `rx` and `rz`.
    this.rx(Math.PI / 2, q);
    this.rz(theta, q);
    this.rx(-Math.PI / 2, q);
    return this;
  }


  // Applies a z gate to the given qubit.
  QuantumCircuit z(int q) {
    // This gate is constructed from `rz`.
    this.rz(Math.PI, q);
    return this;
  }


  // Applies a y gate to the given qubit.
  QuantumCircuit y(int q) {
    // This gate is constructed from `rz` and `x`.
    this.rz(Math.PI, q);
    this.x(q);
    return this;
  }


  /// Applies a measure gate to the given qubit and bit.
  QuantumCircuit measure(int q, int b) {
    if (q >= this.numQubits) {
      throw new IndexOutOfBoundsException("Index for qubit out of range.");
    }
    if (b >= this.numClbits) {
      throw new IndexOutOfBoundsException("Index for output bit out of range.");
    }
    this.data.add(Arrays.asList("m", new Integer(q), new Integer(b)));
    return this;
  }
}


class Simulator {
  double r2 = 0.70710678118;

  //double[][] superpose(int[] x, int[] y) {
  //  double[][] sup = new double[2][2];
  //  sup[0][0] = (x[0] + y[0]) * r2;
  //  sup[0][1] = (x[1] + y[1]) * r2;
  //  sup[1][0] = (x[0] - y[0]) * r2;
  //  sup[1][1] = (x[1] - y[1]) * r2;
  //  return sup;
  //};

  List<List<Double>> superpose(List<Double> x, List<Double> y) {
    List<List<Double>> sup = new ArrayList(Arrays.asList( //<>//
      new ArrayList(Arrays.asList(
        (((Double)x.get(0).doubleValue() + ((Double)y.get(0).doubleValue())) * r2),
        (((Double)x.get(1).doubleValue() + ((Double)y.get(1).doubleValue())) * r2)
      )),
      new ArrayList(Arrays.asList(
        (((Double)x.get(0).doubleValue() - ((Double)y.get(0).doubleValue())) * r2),
        (((Double)x.get(1).doubleValue() - ((Double)y.get(1).doubleValue())) * r2)
      ))
    ));
    return sup; //<>//
  };

  List<List<Double>> turn(List<Double> x, List<Double> y, double theta) {
    List<List<Double>> trn = new ArrayList(Arrays.asList(
      new ArrayList(Arrays.asList(
        ((Double)x.get(0).doubleValue() * Math.cos(theta / 2) + //<>//
        ((Double)y.get(1).doubleValue()) * Math.sin(theta / 2)),
        ((Double)x.get(1).doubleValue() * Math.cos(theta / 2) -
        ((Double)y.get(1).doubleValue()) * Math.sin(theta / 2))
      )),
      new ArrayList(Arrays.asList(
        ((Double)y.get(0).doubleValue() * Math.cos(theta / 2) +
        ((Double)x.get(1).doubleValue()) * Math.sin(theta / 2)),
        ((Double)y.get(1).doubleValue()  * Math.cos(theta / 2) -
        ((Double)x.get(0).doubleValue()) * Math.sin(theta / 2))
      ))
    ));
    return trn;
  };


//  List<List<Double>> simulate(QuantumCircuit qc, int shots, String get) {
  List<List> simulate(QuantumCircuit qc, int shots, String get) {
    List<List> k = new ArrayList();
    for (int j = 0; j < Math.pow(2, qc.numQubits); j++) {
      k.add(Arrays.asList(0.0d, 0.0d));
    }
    k.set(0, Arrays.asList(1.0d, 0.0d));
 //<>//
    Map outputMap = new HashMap();
    for (int idx = 0; idx < qc.data.size(); idx++) {
      List gate = qc.data.get(idx);
      if (gate.get(0).equals("m")) {
        outputMap.put(gate.get(2), gate.get(1));
      }
      else if (gate.get(0).equals("x") || gate.get(0).equals("h") || gate.get(0).equals("rx")) {
        int j = ((Integer)gate.get(gate.size() - 1)).intValue();

        for (int i0 = 0; i0 < Math.pow(2, j); i0++) {
          for (int i1 = 0; i1 < Math.pow(2, qc.numQubits - j - 1); i1++) {
            int b0 = (int)(i0 + Math.pow(2, (j + 1)) * i1);
            int b1 = (int)(b0 + Math.pow(2, j));
            if (gate.get(0).equals("x")) {
              println("b0: " + b0);
              println("b1: " + b1);
              List temp0 = k.get(b0);
              List temp1 = k.get(b1);
              k.set(b0, temp1);
              k.set(b1, temp0);
            }
            else if (gate.get(0).equals("h")) { //<>//
              List<List<Double>> sup = this.superpose(k.get(b0), k.get(b1));
              k.set(b0, sup.get(0));
              k.set(b1, sup.get(1));
            }
            else {
              double theta = ((Double)gate.get(1)).doubleValue();
              List<List<Double>> trn = this.turn(k.get(b0), k.get(b1), theta);
              k.set(b0, trn.get(0));
              k.set(b1, trn.get(1));
            }
          }
        }
      }

      else if (gate.get(0).equals("cx")) {
        int s = ((Integer)gate.get(1)).intValue();
        int t = ((Integer)gate.get(2)).intValue();
        int l = Math.min(s, t);
        int h = Math.max(s, t);
        for (int i0 = 0; i0 < Math.pow(2, l); i0++) {
          for (int i1 = 0; i1 < Math.pow(2, (h - l - 1)); i1++) {
            for (int i2 = 0; i2 < Math.pow(2, (qc.numQubits - h - 1)); i2++) {
              int b0 = (int)(i0 + Math.pow(2, l + 1) * i1 + Math.pow(2, h + 1) * i2 + Math.pow(2, s));
              int b1 = (int)(b0 + Math.pow(2, t));
              List<Double> tmp0 = k.get(b0);
              List<Double> tmp1 = k.get(b1);
              k.set(b0, tmp1);
              k.set(b1, tmp0);
            }
          }
        }
      }
    }
    if (get.equals("statevector")) {
      return k;
    }
    else {
      List<Boolean> m = new ArrayList();
      for (int idx = 0; idx < qc.numQubits; idx++) {
        m.add(false);
      }
      for (int i = 0; i < qc.data.size(); i++) {
        List gate = qc.data.get(i);
        for (int j = 0; j < qc.numQubits; j++) {
          //int j = ((Integer)gate.get(gate.size() - 1)).intValue();
          if (((((Integer)gate.get(gate.size() - 1)).intValue() == j) && m.get(j))) {
            println("Incorrect or missing measure command.");
          }
          m.set(j, gate.get(0).equals("m") && 
            ((Integer)gate.get(1)).intValue() == j && 
            ((Integer)gate.get(2)).intValue() == j);
        }
      }

      List<Double> probs = new ArrayList();
      for (int i = 0; i < k.size(); i++) {
        probs.add((Math.pow(((Double)k.get(i).get(0)).doubleValue(), 2.0d) + 
          Math.pow(((Double)k.get(i).get(1)).doubleValue(), 2.0d)));
      }
      if (get.equals("counts")) {
        /*
        var me = [];
        for (var idx = 0; idx < shots; idx++) {
          var cumu = 0.0;
          var un = true;
          var r = Math.random();
          for (var j = 0; j < probs.length; j++) {
            var p = probs[j];
            cumu += p;
            if (r < cumu && un) {
              var bitStr = j.toString(2);
              var padStr = Math.pow(10, qc.numQubits - bitStr.length).toString().substr(1, qc.numQubits);
              var rawOut = padStr + bitStr;
              var outList = [];
              for (var i = 0; i < qc.numClbits; i++) {
                outList.push('0');
              }
              for (var bit in outputMap) {
                outList[qc.numClbits - 1 - bit] =
                  rawOut[qc.numQubits - 1 - outputMap[bit]];
              }
              var out = outList.join("");
              me.push(out);
              un = false;
            }
          }
        }
        var counts = {};
        for (var meIdx = 0; meIdx < me.length; meIdx++) {
          var out = me[meIdx];
          if (counts.hasOwnProperty(out)) {
            counts[out] += 1;
          } else {
            counts[out] = 1;
          }
        }
        return counts;
        */
      }
    }
    return null;

  }
}


QuantumCircuit psiMinus = new QuantumCircuit(2, 2);
psiMinus.h(0);
psiMinus.h(1);
//psiMinus.x(1);
//psiMinus.cx(0, 1);
//psiMinus.z(1);
psiMinus.measure(0, 0);
psiMinus.measure(1, 1);
List<List> psiMinusStatevector =
  new Simulator().simulate(psiMinus, 0, "statevector");
println("psiMinusStatevector: " + psiMinusStatevector);

List<List> psiMinusCounts =
  new Simulator().simulate(psiMinus, 0, "counts");
println("psiMinusCounts: " + psiMinusCounts);




/*
QuantumCircuit qc = new QuantumCircuit(3, 3);
qc.x(0);
qc.rx(Math.PI, 1);
qc.x(1);
qc.h(2);
qc.cx(0, 1);
qc.z(1);
qc.rz(Math.PI / 2, 1);
qc.ry(Math.PI / 4, 1);
qc.measure(0, 0);
qc.measure(1, 1);
qc.measure(2, 2);
*/

/*
qc.x(0);
qc.rx(Math.PI, 1);
qc.h(1);
qc.cx(2, 1);
qc.rz(Math.PI/2, 1);
qc.ry(Math.PI/4, 1);
qc.measure(1, 1);
*/

//println(psiMinus.numQubits);
//println(psiMinus.numClbits);
//println(psiMinus.data);


//Simulator simulator = new Simulator();
//println(simulator.simulate(qc, 1024, "statevector"));

/*
int[] a = {1, 2};
int[] b = {1, 2};
println(simulator.superpose(a, b).toString());
println(simulator.turn(a, b, Math.PI).toString());
*/
