const accentColor = [50,100,150];
const htmlCanvasId = "canvasPlane";
/* DOM interaction */

htmlVisuals = (function() {
  return {
    updateGroupLabel: function(group, level) {
      groupLabel = document.getElementById("groupLabel")
	    katex.render(group.tex, groupLabel, {throwOnError: false});
    },
    updateHovering: function(moebious) {
      var popup = document.getElementById("hover-popup");
      if(moebious == null) {
        popup.style.display = "none"
        return;
      }
      katex.render(moebious.tex(), popup, {throwOnError: false});
      popup.style.display = "block"
    }
  }
}());

document.addEventListener('mousemove', function(ev) {
  mouseX=ev.pageX;
  mouseY=ev.pageY;
  popup = document.getElementById("hover-popup")
  popup.style.left = (mouseX+10)+"px"
  popup.style.bottom = (window.innerHeight+10-mouseY)+"px"
});

/* Canvas and drawing */

halfplane = (function() {
  var options = {
    mapping: {
      scale: 200,
      origin: [200, 350],
      changed: false,
    },
    hoverInfo: true,
    style: {
      domain: {
        fill: [...accentColor, 100],
        stroke: [...accentColor, 100],
        strokeWeight: 1,
      },
      hover: {
        fill: null,
        stroke: [...accentColor, 150],
        strokeWeight: 3,
      }
    }
  }

  var canvas = null
  var sketch = null

  var domainPrerender = null

  var cosetRepr = {
    group: math.congruenceSubgroups.Gamma_0,
    level: 5,
    matrices: math.congruenceSubgroups.Gamma_0.cosetReprs(5),
    changeLevel: function(level) {
      level = parseInt(level)
      this.matrices = this.group.cosetReprs(level)
      this.level = level
      this.changed = true
      htmlVisuals.updateGroupLabel(this.group, this.level)
    },
    changeGroup: function(group) {
      this.matrices = group.cosetReprs(this.level)
      this.group = group
      this.changed = true
      htmlVisuals.updateGroupLabel(this.group, this.level)
    },
    changed: true
  }

  function prerenderDomain() {
    domainPrerender.clear()
    useStyle(domainPrerender, options.style.domain)
    for (var m of cosetRepr.matrices) {
      drawDomain(domainPrerender, math.congruenceSubgroups.fundamentalDomains[0], m)
    }
  }

  function setup(p) {
    var renderer = p.createCanvas(400, 400)
    domainPrerender = p.createGraphics(400, 400)
    renderer.id(htmlCanvasId)
    p.textFont("Computer Modern")
  }

  function draw(p) {
    canvas.background(255)
    if (cosetRepr.changed || options.mapping.changed) {
      cosetRepr.changed = false
      options.mapping.changed = false
      prerenderDomain()
    }
    canvas.image(domainPrerender, 0, 0)
    drawHovered(canvas)
    drawAxes(canvas)
  }

  function gcd(a, b) {
    if (!b) {
      return a;
    }

    return gcd(b, a % b);
  }

  function drawHovered(g) {

    var p = invMapping(g.mouseX, g.mouseY)
    var M = math.congruenceSubgroups.fundamentalDomains[0].findMoebiousToDomain(p)
    if (p.im >= 0 && M != null && (0 <= g.mouseX && g.mouseX <= g.width && 0 <= g.mouseY && g.mouseY <= g.height)) {
      var ind = cosetRepr.group.cosetReprIndexIn(cosetRepr.level, cosetRepr.matrices, M.inv())
      if (ind != -1) {
        var M2 = cosetRepr.matrices[ind]

        var m = mapping(M2.transform(M.transform(p)))
        g.strokeWeight(5)
        g.point(m[0], m[1])
      }

      var M2 = M.inv()
      useStyle(g, options.style.hover)
      drawDomain(g, math.congruenceSubgroups.fundamentalDomains[0], M2)

      /* Annotate cusp */
      if (M2.m[2] != 0) {
        canvas.strokeWeight(0)
        canvas.fill(0)
        canvas.stroke(0)
        var p = M2.m[0]
        var q = M2.m[2]
        if (p / q != 0 && p / q != 0.5) {

          var d = gcd(q, p)
          if (q / d < 0) {
            q = -q
            p = -p
          }
          annotateX(canvas, p / d, q / d)
        }
      }

      htmlVisuals.updateHovering(M2)
    } else {
      htmlVisuals.updateHovering(null)
    }
  }



  /* Applies fill, stroke and weight */
  function useStyle(graphics, style) {
    if (style.fill == null)
      graphics.noFill()
    else
      graphics.fill(...style.fill)
    if (style.stroke == null || style.strokeWeight == 0 || style.strokeWeight == null)
      graphics.noStroke()
    else {
      graphics.strokeWeight(style.strokeWeight)
      graphics.stroke(...style.stroke)
    }
  }

  function drawDomain(g, domain, m) {
    drawIterator(g, traceDomain(domain, m))
  }

  function* traceDomain(domain, m) {
    var corners = domain.corners
    var i = 1
    for (; i < corners.length; i++)
      yield* mathPainter.traceHyperbolicLine(m.transform(corners[i - 1]), m.transform(corners[i]))
    yield* mathPainter.traceHyperbolicLine(m.transform(corners[i - 1]), m.transform(corners[0]))
  }

  function mapping(cnumber) {
    const s = options.mapping.scale
    const o = options.mapping.origin
    return [cnumber.re * s + o[0], -cnumber.im * s + o[1]]
  }

  function invMapping(pX, pY) {
    const s = options.mapping.scale
    const o = options.mapping.origin
    return new math.Complex((pX - o[0]) / s, (-pY + o[1]) / s)
  }

  /* Iterator through cplx values */
  function drawIterator(g, iter, close = false) {
    g.beginShape()
    var prev = null
    for (var x of iter) {
      if (x != math.oo) {
        var p = mapping(x)
        if (prev == math.oo) {
          g.vertex(p[0], -10)
        }

        g.vertex(p[0], p[1])

      } else if (prev != math.oo && prev != null)
        g.vertex(mapping(prev)[0], -10)
      prev = x
    }
    if (close)
      g.endShape(g.CLOSE)
    else
      g.endShape()
  }

  function annotateX(g, p, q, size = 10) {
    var m = mapping(new math.Complex(1.0 * p / q))
    p = "" + p
    q = "" + q
    g.textSize(size)
    g.textAlign(g.CENTER, g.TOP)
    var w = Math.max(g.textWidth(p), g.textWidth(q))

    g.text(p, m[0], m[1] + 3)
    g.text(q, m[0], m[1] + 3 + size)
    g.strokeWeight(size / 10.0)
    g.line(m[0] - w / 2, m[1] + 2 + size, m[0] + w / 2, m[1] + 2 + size)
    g.line(m[0], m[1], m[0], m[1] + 2)
  }

  function drawAxes(g) {
    const orig = options.mapping.origin
    const width = g.width
    const height = g.height

    g.strokeWeight(1.3);
    g.stroke(0);
    g.fill(0)
    g.line(orig[0], 3, orig[0], height)
    g.line(0, orig[1], width - 3, orig[1])
    g.strokeWeight(0)
    g.triangle(orig[0], 0, orig[0] - 3, 5, orig[0] + 3, 5)
    g.triangle(width, orig[1], width - 5, orig[1] + 3, width - 5, orig[1] - 3)

    g.textStyle(g.NORMAL)
    g.textAlign(g.RIGHT, g.TOP)
    g.textSize(10)
    g.text("Re", width - 3, orig[1] + 3)
    g.text("Im", orig[0] - 3, 3)
    g.text("0", orig[0] - 3, orig[1] + 3)

    annotateX(g, 1, 2)
  }

  sketch = function(p) {
    var dragFromCanvas = false;
    p.setup = function() {
      setup(p)
    };
    p.draw = function() {
      draw(p)
    };
    p.mouseWheel = function(event) {
        if(event.srcElement.id != htmlCanvasId) {
          return;
        }
        var d = event.delta;
        options.mapping.scale = options.mapping.scale*(Math.exp(-0.001*d));
        options.mapping.changed = true;
    };
    p.mouseDragged = function(event) {
      if(!dragFromCanvas) {
        return;
      }
      var d = event.movementX;
      options.mapping.origin[0] = options.mapping.origin[0]+d;
      options.mapping.changed = true;
    };
    p.mousePressed = function(event) {
      dragFromCanvas = event.srcElement.id == htmlCanvasId;
    };
  }
  canvas = new p5(sketch, "drawing")

  return {
    options: options,
    changeGroup: function(group) {
      cosetRepr.changeGroup(group)
    },
    changeLevel: function(level) {
      cosetRepr.changeLevel(level)
    },
    sketch: sketch,
    p5: canvas
  }
}());







//
