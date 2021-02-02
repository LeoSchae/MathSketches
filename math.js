
const math = (function() {

  const InfinityType = "Infinity"
  const ComplexType = "Complex"
  const MoebiousType = "Moebious"

  const oo = {
    mtype: InfinityType
  }

  class Complex {
    constructor(re,im=0) {
      this.mtype = ComplexType
      this.re = re
      this.im = im
    }
    add(val) {
      switch(typeof(val)){
        case 'number':
          return new Complex(this.re+val,this.im)
        default:
          return new Complex(this.re+val.re,this.im+val.im)
      }
    }
    sub(val) {
      switch(typeof(val)){
        case 'number':
          return new Complex(this.re-val,this.im)
        default:
          return new Complex(this.re-val.re,this.im-val.im)
      }
    }
    mul(val) {
      switch(typeof(val)){
        case 'number':
          return new Complex(val*this.re,val*this.im)
        default:
          return new Complex(this.re*val.re-this.im*val.im,this.im*val.re+this.re*val.im)
      }
    }
    div(val) {
      switch(typeof(val)){
        case 'number':
          return new Complex(this.re/val,this.im/val)
        default:
          return this.mul(val.inv())
      }
    }
    abs2() {
      return this.re*this.re+this.im*this.im
    }
    abs() {
      return Math.sqrt(this.abs2())
    }
    inv() {
      const n = this.abs2()
      return new Complex(this.re/n,-this.im/n)
    }
    polarAngle() {
      if (this.re == 0 && this.im == 0) {
        return 0;
      }
      if(this.re == 0)
        return Math.PI/2
      if (this.re >= 0) {
        const phi = Math.atan(1.0 * this.im / this.re);
        if (phi < 0)
          return 2 * Math.PI + phi;
        return phi;
      }
      return Math.atan(1.0 * this.im / this.re) + Math.PI;
    }
  }

  /* Only works for ad-bc = 1. Otherwise the inverse breaks */
  class Moebious {
    constructor(a, b, c, d) {
      this.mtype = MoebiousType
      this.m = [a, b, c, d]
    }
    mul(val) {
      const m1 = this.m
      const m2 = val.m
      return new Moebious(m1[0]*m2[0]+m1[1]*m2[2],m1[0]*m2[1]+m1[1]*m2[3],m1[2]*m2[0]+m1[3]*m2[2],m1[2]*m2[1]+m1[3]*m2[3])
    }
    transform(val) {
      switch(typeof(val)) {
        case 'number':
          val = new Complex(val)
        default:
          if (val == oo) {
            if (this.m[2] == 0)
              return oo
            return new Complex(this.m[0] / this.m[2])
          }
          const q = val.mul(this.m[2]).add(this.m[3])
          if (q.re == 0 && q.im == 0)
            return oo
          return val.mul(this.m[0]).add(this.m[1]).div(q)
      }
    }
    inv() {
      var m = this.m
      return new Moebious(m[3],-m[1],-m[2],m[0])
    }
    tex() {
      var m = this.m
      return "\\begin{pmatrix}"+m[0]+"&"+m[1]+"\\\\"+m[2]+"&"+m[3]+"\\end{pmatrix}"
    }
  }



  const congruenceSubgroups = (function(){
    class CongruenceSubgroup {
      /**
       * tex: the name given in TeX
       * indicatorFunction: a function (level,value) -> true/false whether value is in the group
       */
      constructor(tex,indicatorFunction) {
        this.tex = tex
        this.isInGroup = indicatorFunction
      }
      cosetReprIndexIn(level, list, value) {
        var xinv = value.inv();

        for (var i = 0, l=list.length; i < l; i++) {
          if (this.isInGroup(level,list[i].mul(xinv))) {
            return i;
          }
        }
        return -1;
      }
      cosetReprIn(level,list,value) {
        var i = this.cosetReprIndexIn(level,list,value)
        return i==-1 ? null : list[i]
      }
      cosetReprs(level) {
        const generators = [new Moebious(0, -1, 1, 0), new Moebious(1, 1, 0, 1), new Moebious(1, -1, 0, 1)]
        const group = this

        var completeList = [new Moebious(1, 0, 0, 1)]

        var checks = [] // 2 iterations ago
        var seeds = []  // 1 iteration  ago
        var added = [new Moebious(1, 0, 0, 1)]  // current iteration
        var it = 0
        while(added.length > 0) {
          checks = seeds
          seeds = added
          added = []

          for(var s of seeds) {
            for(var g of generators) {
              var x = s.mul(g) // step in generator direction
              if(group.cosetReprIndexIn(level,checks,x) == -1 && group.cosetReprIndexIn(level,added,x) == -1 && group.cosetReprIndexIn(level,seeds,x) == -1) {

                completeList.push(x);
                added.push(x);
              }
            }
          }

        }
        return completeList;
      }
    }

    function modN(x,N) {
      return ((x % N)+N)%N
    }

    function gamma_0_indicator(level,value) {
      const m = value.m
      return modN(m[2],level) == 0
    }

    function gamma_1_indicator(level,value) {
      const m = value.m
      return modN(m[2],level) == 0 && (modN(m[0],level)==1 || modN(m[0],level)==level-1)
    }

    function gamma_indicator(level,value) {
      const m = value.m
      return modN(m[2],level) == 0 && modN(m[1],level) == 0 && (modN(m[0],level)==1 || modN(m[0],level)==level-1)
    }

    function exp_i(x) {
      return new Complex(Math.cos(x), Math.sin(x))
    }

    domain1 = {
      corners: [exp_i(Math.PI / 3), oo, exp_i(2 * Math.PI / 3)],
      findMoebiousToDomain: function(p, maxIter = 1000) {
        var result = new Moebious(1, 0, 0, 1)
        var current = p;

        for (var i = 0; i < maxIter; i++) {
          if (current == oo)
            return result

          var n = Math.floor(current.re + 0.5)
          result = new Moebious(1, -n, 0, 1).mul(result)
          current = result.transform(p)

          if (current.abs2() >= 1) {
            return result;
          }

          result = new Moebious(0, 1, -1, 0).mul(result)
          current = result.transform(p)
        }
        console.log("Max iterations in 'findMoebiousToFund' reached")
        return null
      }
    }
    domain2 = {
      corners: [new Complex(0, 1), oo, new Complex(1, 1), exp_i(Math.PI / 3)],
      findMoebiousToDomain: function(p, maxIter=1000) {
        var result = domain1.findMoebiousToFund1(p,maxIter)
        if (result == null)
          return result
        var current = result.transform(p)

        // If -0.5 <= Re(z) < 0 shift by z+1 to the right domain
        if (current.re < 0)
          return new Moebious(1, 1, 0, 1).mul(result)
        return result
      }
    }

    return {
      Gamma_0: new CongruenceSubgroup("\\Gamma_0",gamma_0_indicator),
      Gamma_1: new CongruenceSubgroup("\\Gamma_1",gamma_1_indicator),
      Gamma:   new CongruenceSubgroup("\\Gamma",gamma_indicator),
      fundamentalDomains: [domain1,domain2]
    }
  })()


  /* Parses all objects to mtypes if possible */
  function parse(object) {
    if(Array.isArray(object)) {
      const l = object.length
      for(var i = 0; i<l; i++) {
        object[i] = parse(object[i])
      }
      return object
    } else {
      if("mtype" in object)
        return asMType(object)
    }
  }

  function asMType(object) {
    switch (object.mtype) {
      case ComplexType:
        return new Complex(object.re, object.im)
      case MoebiousType:
        const m = object.m
        return new Moebious(m[0], m[1], m[2], m[3])
      case InfinityType:
        return oo;
      default:
        console.log("Failed to parse MType")
        return undefined
    }
  }

  return {
    oo: oo,
    Complex: Complex,
    Moebious: Moebious,
    parse: parse,
    congruenceSubgroups: congruenceSubgroups
  }
}());

/* Functions for transforming shapes to lines or outlines */
mathPainter = (function(math) {
  const {oo,Complex,Moebious} = math

  function exp_i(x) {
    return new Complex(Math.cos(x), Math.sin(x))
  }

  /**
   * Draw circle arc in normal metric
   */
  function* traceCircleArc(radius, center, startAngle, stopAngle, segments = 15) {
    var s = (stopAngle - startAngle) / segments
    radius = new Complex(radius)
    for (var i = 0; i <= segments; i++)
      yield center.add(radius.mul(exp_i(startAngle + i * s)))
  }

  /**
   * Draw a hyperbolic line between points
   * Segment amount is only used for arcs, not lines.
   */
  function* traceHyperbolicLine(from, to, segments = 15) {
    if (from == oo || to == oo || Math.abs(from.re - to.re) < 0.0000000001) {
      yield from
      yield to
      return
    }
    var c = new Complex((from.abs2() - to.abs2()) / 2 / (from.re - to.re))
    var cf = from.sub(c)
    var ct = to.sub(c)

    yield* traceCircleArc(cf.abs(), c, cf.polarAngle(), ct.polarAngle(), segments)
  }

  return {
    traceCircleArc: traceCircleArc,
    traceHyperbolicLine: traceHyperbolicLine
  }
})(math);





//
