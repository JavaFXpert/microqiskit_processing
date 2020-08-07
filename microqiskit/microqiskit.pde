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
    this.data.add(Arrays.asList("x"));
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
  double[][] superpose(int[] x, int[] y) {
    double[][] sup = new double[2][2];
    sup[0][0] = (x[0] + y[0]) * r2;
    sup[0][1] = (x[1] + y[1]) * r2;
    sup[1][0] = (x[0] - y[0]) * r2;
    sup[1][1] = (x[1] - y[1]) * r2;
    return sup;
  };

  double[][] turn(int[] x, int[] y, double theta) {
    double[][] trn = new double[2][2];
    trn[0][0] = x[0] * Math.cos(theta / 2) + y[1] * Math.sin(theta / 2);
    trn[0][1] = x[1] * Math.cos(theta / 2) - y[0] * Math.sin(theta / 2);
    trn[1][0] = y[0] * Math.cos(theta / 2) + x[1] * Math.sin(theta / 2);
    trn[1][1] = y[1] * Math.cos(theta / 2) - x[0] * Math.sin(theta / 2);
    return trn;
  };
  
  
  List simulate(QuantumCircuit qc, int shots, String get) {
    List<List> k = new ArrayList();
    for (int j = 0; j < Math.pow(2, qc.numQubits); j++) {
      k.add(Arrays.asList(0.0, 0.0));
    }
    k.set(0, Arrays.asList(1.0, 0.0));
    
    Map outputMap = new HashMap();
    for (int idx = 0; idx < qc.data.size(); idx++) {
      List gate = qc.data.get(idx);
      if (gate.get(0).equals("m")) {
        outputMap.set(gate.get(2), gate.get(1);
      }  
      else if (gate.get(0).equals("x") || gate.get(0).equals("h") || gate.get(0).equals("rx")) {
      var j = gate.slice(-1)[0];
      for (var i0 = 0; i0 < Math.pow(2, j); i0++) {
        for (var i1 = 0; i1 < Math.pow(2, qc.numQubits - j - 1); i1++) {
          var b0 = i0 + Math.pow(2, (j + 1)) * i1;
          var b1 = b0 + Math.pow(2, j);
          if (gate[0] == 'x') {
            var temp0 = k[b0];
            var temp1 = k[b1];
            k[b0] = temp1;
            k[b1] = temp0;
          } else if (gate[0] == 'h') {
            var sup = superpose(k[b0], k[b1]);
            k[b0] = sup[0];
            k[b1] = sup[1];
          } else {
            var theta = gate[1];
            var trn = turn(k[b0], k[b1], theta);
            k[b0] = trn[0];
            k[b1] = trn[1];
          }
        }
      }
    }
    
/*
  var k = [];
  for (j = 0; j < Math.pow(2, qc.numQubits); j++) {
    k.push([0, 0]);
  }
  k[0] = [1.0, 0.0];
  var outputMap = {};
  
  for (var idx = 0; idx < qc.data.length; idx++) {
    var gate = qc.data[idx];
    if (gate[0] == 'm') {
      outputMap[gate[2]] = gate[1];
    } else if (gate[0] == "x" || gate[0] == "h" || gate[0] == "rx") {
      var j = gate.slice(-1)[0];
      for (var i0 = 0; i0 < Math.pow(2, j); i0++) {
        for (var i1 = 0; i1 < Math.pow(2, qc.numQubits - j - 1); i1++) {
          var b0 = i0 + Math.pow(2, (j + 1)) * i1;
          var b1 = b0 + Math.pow(2, j);
          if (gate[0] == 'x') {
            var temp0 = k[b0];
            var temp1 = k[b1];
            k[b0] = temp1;
            k[b1] = temp0;
          } else if (gate[0] == 'h') {
            var sup = superpose(k[b0], k[b1]);
            k[b0] = sup[0];
            k[b1] = sup[1];
          } else {
            var theta = gate[1];
            var trn = turn(k[b0], k[b1], theta);
            k[b0] = trn[0];
            k[b1] = trn[1];
          }
        }
      }
    }
    else if (gate[0] == 'cx') {
      var s = gate[1];
      var t = gate[2];
      var l = Math.min(s, t);
      var h = Math.max(s, t);
      for (var i0 = 0; i0 < Math.pow(2, l); i0++) {
        for (var i1 = 0; i1 < Math.pow(2, (h - l - 1)); i1++) {
          for (var i2 = 0; i2 < Math.pow(2, (qc.numQubits - h - 1)); i2++) {
            var b0 = i0 + Math.pow(2, l + 1) * i1 + Math.pow(2, h + 1) * i2 + Math.pow(2, s);
            var b1 = b0 + Math.pow(2, t);
            var tmp0 = k[b0];
            var tmp1 = k[b1];
            k[b0] = tmp1;
            k[b1] = tmp0;
          }
        }
      }
    }
  }
  if (get == 'statevector') {
    return k;
  }
  else {
    var m = [];
    for (var idx = 0; idx < qc.numQubits; idx++) {
      m.push(false);
    }
    for (var i = 0; i < qc.data.length; i++) {
      var gate = qc.data[i];
      for (var j = 0; j < qc.numQubits; j++) {
        if (((gate.slice(-1)[0] == j) && m[j])) {
          throw ('Incorrect or missing measure command.');
        }
        m[j] = (gate[0] == 'm' && gate[1] == j && gate[2] == j);
      }
    }
    var probs = [];
    for (var i = 0; i < k.length; i++) {
      probs.push((Math.pow(k[i][0], 2) + Math.pow(k[i][1], 2)));
    }
    if (get == 'counts' || get == 'memory') {
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
      if (get == 'memory') {
        return m;
      } else {
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
      }
    }
  }
*/
    
    return k;
  }
}



QuantumCircuit qc = new QuantumCircuit(3, 3);
qc.x(0);
qc.rx(Math.PI, 1);    
qc.h(1);
qc.cx(2, 1);
qc.rz(Math.PI/2, 1);
qc.ry(Math.PI/4, 1);
qc.measure(1, 1);
println(qc.numQubits);
println(qc.numClbits);
println(qc.data);

Simulator simulator = new Simulator();
println(simulator.simulate(qc, 1024, "statevector"));

int[] a = {1, 2};
int[] b = {1, 2};
println(simulator.superpose(a, b).toString());
println(simulator.turn(a, b, Math.PI).toString());
