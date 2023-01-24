//'use strict';

//module.exports = {
 //   compactNode: compactNode,
 //   compactGraph: compactGraph
//};

function findNextEnd(prev, v, vertices, ends, vertexCoords, edgeData, trackIncoming, options) {
    var weight = vertices[prev][v],
        reverseWeight = vertices[v][prev],
        coordinates = [],
        path = [],
        reducedEdge = options.edgeDataSeed;
        
    if (options.edgeDataReduceFn) {
        reducedEdge = options.edgeDataReduceFn(reducedEdge, edgeData[v][prev]);
    }

    while (!ends[v]) {
        var edges = vertices[v];

        if (!edges) { break; }

        var next = Object.keys(edges).filter(function notPrevious(k) { return k !== prev; })[0];
        weight += edges[next];

        if (trackIncoming) {
            reverseWeight += vertices[next][v];

            if (path.indexOf(v) >= 0) {
                ends[v] = vertices[v];
                break;
            }
            path.push(v);
        }

        if (options.edgeDataReduceFn) {
            reducedEdge = options.edgeDataReduceFn(reducedEdge, edgeData[v][next]);
        }

        coordinates.push(vertexCoords[v]);
        prev = v;
        v = next;
    }

    return {
        vertex: v,
        weight: weight,
        reverseWeight: reverseWeight,
        coordinates: coordinates,
        reducedEdge: reducedEdge
    };
}

function compactNode(k, vertices, ends, vertexCoords, edgeData, trackIncoming, options) {
    options = options || {};
    var neighbors = vertices[k];
    return Object.keys(neighbors).reduce(function compactEdge(result, j) {
        var neighbor = findNextEnd(k, j, vertices, ends, vertexCoords, edgeData, trackIncoming, options);
        var weight = neighbor.weight;
        var reverseWeight = neighbor.reverseWeight;
        if (neighbor.vertex !== k) {
            if (!result.edges[neighbor.vertex] || result.edges[neighbor.vertex] > weight) {
                result.edges[neighbor.vertex] = weight;
                result.coordinates[neighbor.vertex] = [vertexCoords[k]].concat(neighbor.coordinates);
                result.reducedEdges[neighbor.vertex] = neighbor.reducedEdge;
            }
            if (trackIncoming && 
                !isNaN(reverseWeight) && (!result.incomingEdges[neighbor.vertex] || result.incomingEdges[neighbor.vertex] > reverseWeight)) {
                result.incomingEdges[neighbor.vertex] = reverseWeight;
                var coordinates = [vertexCoords[k]].concat(neighbor.coordinates);
                coordinates.reverse();
                result.incomingCoordinates[neighbor.vertex] = coordinates;
            }
        }
        return result;
    }, {edges: {}, incomingEdges: {}, coordinates: {}, incomingCoordinates: {}, reducedEdges: {}});
}

function compactGraph(vertices, vertexCoords, edgeData, options) {
    options = options || {};
    var progress = options.progress;
    var ends = Object.keys(vertices).reduce(function findEnds(es, k, i, vs) {
        var vertex = vertices[k];
        var edges = Object.keys(vertex);
        var numberEdges = edges.length;
        var remove;

        if (numberEdges === 1) {
            var other = vertices[edges[0]];
            remove = !other[k];
        } else if (numberEdges === 2) {
            remove = edges.filter(function(n) {
                return vertices[n][k];
            }).length === numberEdges;
        } else {
            remove = false;
        }

        if (!remove) {
            es[k] = vertex;
        }

        if (i % 1000 === 0 && progress) {
            progress('compact:ends', i, vs.length);
        }

        return es;
    }, {});

    return Object.keys(ends).reduce(function compactEnd(result, k, i, es) {
        var compacted = compactNode(k, vertices, ends, vertexCoords, edgeData, false, options);
        result.graph[k] = compacted.edges;
        result.coordinates[k] = compacted.coordinates;

        if (options.edgeDataReduceFn) {
            result.reducedEdges[k] = compacted.reducedEdges;
        }

        if (i % 1000 === 0 && progress) {
            progress('compact:nodes', i, es.length);
        }

        return result;
    }, {graph: {}, coordinates: {}, reducedEdges: {}});
};

















//TinyQueue


var Queue= 
function Queue(data, compare) {
    if ( data === void 0 ) data = [];
    if ( compare === void 0 ) compare = defaultCompare;

    this.data = data;
    this.length = this.data.length;
    this.compare = compare;

    if (this.length > 0) {
        for (var i = (this.length >> 1) - 1; i >= 0; i--) { this._down(i); }
    }
};

Queue.prototype.push = function push (item) {
    this.data.push(item);
    this.length++;
    this._up(this.length - 1);
};

Queue.prototype.pop = function pop () {
    if (this.length === 0) { return undefined; }

    var top = this.data[0];
    var bottom = this.data.pop();
    this.length--;

    if (this.length > 0) {
        this.data[0] = bottom;
        this._down(0);
    }

    return top;
};

Queue.prototype.peek = function peek () {
    return this.data[0];
};

Queue.prototype._up = function _up (pos) {
    var ref = this;
        var data = ref.data;
        var compare = ref.compare;
    var item = data[pos];

    while (pos > 0) {
        var parent = (pos - 1) >> 1;
        var current = data[parent];
        if (compare(item, current) >= 0) { break; }
        data[pos] = current;
        pos = parent;
    }

    data[pos] = item;
};

Queue.prototype._down = function _down (pos) {
    var ref = this;
        var data = ref.data;
        var compare = ref.compare;
    var halfLength = this.length >> 1;
    var item = data[pos];

    while (pos < halfLength) {
        var left = (pos << 1) + 1;
        var best = data[left];
        var right = left + 1;

        if (right < this.length && compare(data[right], best) < 0) {
            left = right;
            best = data[right];
        }
        if (compare(best, item) >= 0) { break; }

        data[pos] = best;
        pos = left;
    }

    data[pos] = item;
};

function defaultCompare(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
}












//dijkstra

//var Queue = require('tinyqueue');

//module.exports =

 var findPath = function findPath(graph, start, end) {
    var costs = {};
    costs[start] = 0;
    var initialState = [0, [start], start];
    var queue = new Queue([initialState], function(a, b) { return a[0] - b[0]; });
	
    var explored = {};

    while (queue.length) {
        var state = queue.pop();
        var cost = state[0];
        var node = state[2];
        if (node === end) {
            return state.slice(0, 2);
        }

        var neighbours = graph[node];
        Object.keys(neighbours).forEach(function(n) {
            var newCost = cost + neighbours[n];
            if (!(n in costs) || newCost < costs[n]) {
                costs[n] = newCost;
                var newState = [newCost, state[1].concat([n]), n];
                queue.push(newState);
            }
        });
    }

    return null;
};












//'use strict';

//var findPath = require('./dijkstra'),
 //   preprocess = require('./preprocessor'),
  //  compactor = require('./compactor'),
   // roundCoord = require('./round-coord');

//module.exports = PathFinder;

 PathFinder = function PathFinder(graph, options) {    
    options = options || {};

    if (!graph.compactedVertices) {
        graph = preprocess(graph, options);
		console.log("no graph .cmpactedvertices");
		console.log("graph with preprocess",graph);
    }

    this._graph = graph;
    this._keyFn = options.keyFn || function(c) {
        return c.join(',');
    };
    this._precision = options.precision || 1e-5;
    this._options = options;

    if (Object.keys(this._graph.compactedVertices).filter(function(k) { return k !== 'edgeData'; }).length === 0) {
        throw new Error('Compacted graph contains no forks (topology has no intersections).');
    }
}

PathFinder.prototype = {
    findPath: function(a, b) {
        var start = this._keyFn(roundCoord(a.geometry.coordinates, this._precision)),
            finish = this._keyFn(roundCoord(b.geometry.coordinates, this._precision));

        // We can't find a path if start or finish isn't in the
        // set of non-compacted vertices
        if (!this._graph.vertices[start] || !this._graph.vertices[finish]) {
            return null;
        }

        var phantomStart = this._createPhantom(start);
        var phantomEnd = this._createPhantom(finish);

        var path = findPath(this._graph.compactedVertices, start, finish);

        if (path) {
            var weight = path[0];
            path = path[1];
            return {
                path: path.reduce(function buildPath(cs, v, i, vs) {
                    if (i > 0) {
                        cs = cs.concat(this._graph.compactedCoordinates[vs[i - 1]][v]);
                    }

                    return cs;
                }.bind(this), []).concat([this._graph.sourceVertices[finish]]),
                weight: weight,
                edgeDatas: this._graph.compactedEdges 
                    ? path.reduce(function buildEdgeData(eds, v, i, vs) {
                        if (i > 0) {
                            eds.push({
                                reducedEdge: this._graph.compactedEdges[vs[i - 1]][v]
                            });
                        }

                        return eds;
                    }.bind(this), [])
                    : undefined
            };
        } else {
            return null;
        }

        this._removePhantom(phantomStart);
        this._removePhantom(phantomEnd);
    },

    serialize: function() {
        return this._graph;
    },

    _createPhantom: function(n) {
        if (this._graph.compactedVertices[n]) return null;

        var phantom = compactNode(n, this._graph.vertices, this._graph.compactedVertices, this._graph.sourceVertices, this._graph.edgeData, true, this._options);
        this._graph.compactedVertices[n] = phantom.edges;
        this._graph.compactedCoordinates[n] = phantom.coordinates;

        if (this._graph.compactedEdges) {
            this._graph.compactedEdges[n] = phantom.reducedEdges;
        }

        Object.keys(phantom.incomingEdges).forEach(function(neighbor) {
            this._graph.compactedVertices[neighbor][n] = phantom.incomingEdges[neighbor];
            this._graph.compactedCoordinates[neighbor][n] = [this._graph.sourceVertices[neighbor]].concat(phantom.incomingCoordinates[neighbor].slice(0, -1));
            if (this._graph.compactedEdges) {
                this._graph.compactedEdges[neighbor][n] = phantom.reducedEdges[neighbor];
            }
        }.bind(this));

        return n;
    },

    _removePhantom: function(n) {
        if (!n) return;

        Object.keys(this._graph.compactedVertices[n]).forEach(function(neighbor) {
            delete this._graph.compactedVertices[neighbor][n];
        }.bind(this));
        Object.keys(this._graph.compactedCoordinates[n]).forEach(function(neighbor) {
            delete this._graph.compactedCoordinates[neighbor][n];
        }.bind(this));
        if (this._graph.compactedEdges) {
            Object.keys(this._graph.compactedEdges[n]).forEach(function(neighbor) {
                delete this._graph.compactedEdges[neighbor][n];
            }.bind(this));
        }

        delete this._graph.compactedVertices[n];
        delete this._graph.compactedCoordinates[n];

        if (this._graph.compactedEdges) {
            delete this._graph.compactedEdges[n];
        }
    }
};









//PREPROCESSOR

//'use strict';

//var topology = require('./topology'),
 //   compactor = require('./compactor'),
  //  distance = require('@turf/distance').default,
  //  roundCoord = require('./round-coord'),
 //   point = require('turf-point');

//module.exports =
 function preprocess(graph, options) {
    options = options || {};
    var weightFn = options.weightFn || function defaultWeightFn(a, b) {
            return distance(point(a), point(b));
        },
        topo;

    if (graph.type === 'FeatureCollection') {
        // Graph is GeoJSON data, create a topology from it
        topo = topology(graph, options);
		//alert("FeatureCollection");
    } 
	else if (graph.edges) {
        // Graph is a preprocessed topology
        topo = graph;
			//alert("graph");
    }

    var graph = topo.edges.reduce(function buildGraph(g, edge, i, es) {
        var a = edge[0],
            b = edge[1],
            props = edge[2],
            w = weightFn(topo.vertices[a], topo.vertices[b], props),
            makeEdgeList = function makeEdgeList(node) {
                if (!g.vertices[node]) {
                    g.vertices[node] = {};
                    if (options.edgeDataReduceFn) {
                        g.edgeData[node] = {};
                    }
                }
            },
            concatEdge = function concatEdge(startNode, endNode, weight) {
                var v = g.vertices[startNode];
                v[endNode] = weight;
                if (options.edgeDataReduceFn) {
                    g.edgeData[startNode][endNode] = options.edgeDataReduceFn(options.edgeDataSeed, props);
                }
            };

        if (w) {
            makeEdgeList(a);
            makeEdgeList(b);
            if (w instanceof Object) {
                if (w.forward) {
                    concatEdge(a, b, w.forward);
                }
                if (w.backward) {
                    concatEdge(b, a, w.backward);
                }
            } else {
                concatEdge(a, b, w);
                concatEdge(b, a, w);
            }
        }

        if (i % 1000 === 0 && options.progress) {
            options.progress('edgeweights', i,es.length);
        }

        return g;
    }, {edgeData: {}, vertices: {}});

    var compact = compactGraph(graph.vertices, topo.vertices, graph.edgeData, options);
console.log("compact", compact);
    return {
        vertices: graph.vertices,
        edgeData: graph.edgeData,
        sourceVertices: topo.vertices,
        compactedVertices: compact.graph,
        compactedCoordinates: compact.coordinates,
        compactedEdges: options.edgeDataReduceFn ? compact.reducedEdges : null
    };
	

}







//roundCoord

//module.exports = 
function roundCoord(c, precision) {
    return [
        Math.round(c[0] / precision) * precision,
        Math.round(c[1] / precision) * precision,
    ];
};





















//'use strict';

//var explode = require('@turf/explode'),
//    roundCoord = require('./round-coord');

//module.exports = topology;

 function geoJsonReduce(geojson, fn, seed) {
    if (geojson.type === 'FeatureCollection') {
        return geojson.features.reduce(function reduceFeatures(a, f) {
            return geoJsonReduce(f, fn, a);
        }, seed);
    } else {
        return fn(seed, geojson);
    }
}

function geoJsonFilterFeatures(geojson, fn) {
    var features = [];
    if (geojson.type === 'FeatureCollection') {
        features = features.concat(geojson.features.filter(fn));
    }

    return {
        type: 'FeatureCollection',
        features: features
    };
}

function isLineString(f) {
    return f.geometry.type === 'LineString';
}

function topology(geojson, options) {
    options = options || {};
    var keyFn = options.keyFn || function defaultKeyFn(c) {
            return c.join(',');
        },
        precision = options.precision || 1e-5;

    var lineStrings = geoJsonFilterFeatures(geojson, isLineString);
    var explodedLineStrings = explode(lineStrings);
    var vertices = explodedLineStrings.features.reduce(function buildTopologyVertices(cs, f, i, fs) {
            var rc = roundCoord(f.geometry.coordinates, precision);
            cs[keyFn(rc)] = f.geometry.coordinates;

            if (i % 1000 === 0 && options.progress) {
                options.progress('topo:vertices', i, fs.length);
            }

            return cs;
        }, {}),
        edges = geoJsonReduce(lineStrings, function buildTopologyEdges(es, f, i, fs) {
            f.geometry.coordinates.forEach(function buildLineStringEdges(c, i, cs) {
                if (i > 0) {
                    var k1 = keyFn(roundCoord(cs[i - 1], precision)),
                        k2 = keyFn(roundCoord(c, precision));
                    es.push([k1, k2, f.properties]);
                }
            });

            if (i % 1000 === 0 && options.progress) {
                options.progress('topo:edges', i, fs.length);
            }

            return es;
        }, []);

    return {
        vertices: vertices,
        edges: edges
    };
}





















//var featureCollection = require('turf-helpers').featureCollection;
//var each = require('turf-meta').coordEach;
//var point = require('turf-helpers').point;

/**
 * Takes a feature or set of features and returns all positions as
 * {@link Point|points}.
 *
 * @name explode
 * @param {(Feature|FeatureCollection)} geojson input features
 * @return {FeatureCollection<point>} points representing the exploded input features
 * @throws {Error} if it encounters an unknown geometry type
 * @example
 * var poly = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Polygon",
 *     "coordinates": [[
 *       [177.434692, -17.77517],
 *       [177.402076, -17.779093],
 *       [177.38079, -17.803937],
 *       [177.40242, -17.826164],
 *       [177.438468, -17.824857],
 *       [177.454948, -17.796746],
 *       [177.434692, -17.77517]
 *     ]]
 *   }
 * };
 *
 * var points = turf.explode(poly);
 *
 * //=poly
 *
 * //=points
 */
//module.exports = 

each= coordEach;
function explode(geojson) {
    var points = [];
    each(geojson, function (coord) {
        points.push(point(coord));
    });
    return featureCollection(points);
};
























//"use strict";
//Object.defineProperty(exports, "__esModule", { value: true });
/**
 * @module helpers
 */
/**
 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
 *
 * @memberof helpers
 * @type {number}
 */
// 
earthRadius = 6371008.8;
/**
 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
 *
 * @memberof helpers
 * @type {Object}
 */
// 
factors = {
    centimeters:  earthRadius * 100,
    centimetres:  earthRadius * 100,
    degrees:  earthRadius / 111325,
    feet:  earthRadius * 3.28084,
    inches:  earthRadius * 39.37,
    kilometers:  earthRadius / 1000,
    kilometres:  earthRadius / 1000,
    meters:  earthRadius,
    metres:  earthRadius,
    miles:  earthRadius / 1609.344,
    millimeters:  earthRadius * 1000,
    millimetres:  earthRadius * 1000,
    nauticalmiles:  earthRadius / 1852,
    radians: 1,
    yards:  earthRadius * 1.0936,
};
/**
 * Units of measurement factors based on 1 meter.
 *
 * @memberof helpers
 * @type {Object}
 */
// 
unitsFactors = {
    centimeters: 100,
    centimetres: 100,
    degrees: 1 / 111325,
    feet: 3.28084,
    inches: 39.37,
    kilometers: 1 / 1000,
    kilometres: 1 / 1000,
    meters: 1,
    metres: 1,
    miles: 1 / 1609.344,
    millimeters: 1000,
    millimetres: 1000,
    nauticalmiles: 1 / 1852,
    radians: 1 /  earthRadius,
    yards: 1.0936133,
};
/**
 * Area of measurement factors based on 1 square meter.
 *
 * @memberof helpers
 * @type {Object}
 */
// 
areaFactors = {
    acres: 0.000247105,
    centimeters: 10000,
    centimetres: 10000,
    feet: 10.763910417,
    hectares: 0.0001,
    inches: 1550.003100006,
    kilometers: 0.000001,
    kilometres: 0.000001,
    meters: 1,
    metres: 1,
    miles: 3.86e-7,
    millimeters: 1000000,
    millimetres: 1000000,
    yards: 1.195990046,
};
/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature} a GeoJSON Feature
 * @example
 * var geometry = {
 *   "type": "Point",
 *   "coordinates": [110, 50]
 * };
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geom, properties, options) {
    if (options === void 0) { options = {}; }
    var feat = { type: "Feature" };
    if (options.id === 0 || options.id) {
        feat.id = options.id;
    }
    if (options.bbox) {
        feat.bbox = options.bbox;
    }
    feat.properties = properties || {};
    feat.geometry = geom;
    return feat;
}
// feature = feature;
/**
 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
 * For GeometryCollection type use `helpers.geometryCollection`
 *
 * @name geometry
 * @param {string} type Geometry Type
 * @param {Array<any>} coordinates Coordinates
 * @param {Object} [options={}] Optional Parameters
 * @returns {Geometry} a GeoJSON Geometry
 * @example
 * var type = "Point";
 * var coordinates = [110, 50];
 * var geometry = turf.geometry(type, coordinates);
 * // => geometry
 */
function geometry(type, coordinates, _options) {
    if (_options === void 0) { _options = {}; }
    switch (type) {
        case "Point":
            return point(coordinates).geometry;
        case "LineString":
            return lineString(coordinates).geometry;
        case "Polygon":
            return polygon(coordinates).geometry;
        case "MultiPoint":
            return multiPoint(coordinates).geometry;
        case "MultiLineString":
            return multiLineString(coordinates).geometry;
        case "MultiPolygon":
            return multiPolygon(coordinates).geometry;
        default:
            throw new Error(type + " is invalid");
    }
}
// geometry = geometry;
/**
 * Creates a {@link Point} {@link Feature} from a Position.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Point>} a Point feature
 * @example
 * var point = turf.point([-75.343, 39.984]);
 *
 * //=point
 */
function point(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    if (!coordinates) {
        throw new Error("coordinates is required");
    }
    if (!Array.isArray(coordinates)) {
        throw new Error("coordinates must be an Array");
    }
    if (coordinates.length < 2) {
        throw new Error("coordinates must be at least 2 numbers long");
    }
    if (!isNumber(coordinates[0]) || !isNumber(coordinates[1])) {
        throw new Error("coordinates must contain numbers");
    }
    var geom = {
        type: "Point",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// point = point;
/**
 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
 *
 * @name points
 * @param {Array<Array<number>>} coordinates an array of Points
 * @param {Object} [properties={}] Translate these properties to each Feature
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Point>} Point Feature
 * @example
 * var points = turf.points([
 *   [-75, 39],
 *   [-80, 45],
 *   [-78, 50]
 * ]);
 *
 * //=points
 */
function points(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return point(coords, properties);
    }), options);
}
// points = points;
/**
 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Polygon>} Polygon Feature
 * @example
 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
 *
 * //=polygon
 */
function polygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    for (var _i = 0, coordinates_1 = coordinates; _i < coordinates_1.length; _i++) {
        var ring = coordinates_1[_i];
        if (ring.length < 4) {
            throw new Error("Each LinearRing of a Polygon must have 4 or more Positions.");
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            // Check if first point of Polygon contains two numbers
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error("First and last Position are not equivalent.");
            }
        }
    }
    var geom = {
        type: "Polygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// polygon = polygon;
/**
 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
 *
 * @name polygons
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
 * @example
 * var polygons = turf.polygons([
 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
 * ]);
 *
 * //=polygons
 */
function polygons(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return polygon(coords, properties);
    }), options);
}
// polygons = polygons;
/**
 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<LineString>} LineString Feature
 * @example
 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
 *
 * //=linestring1
 * //=linestring2
 */
function lineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    if (coordinates.length < 2) {
        throw new Error("coordinates must be an array of two or more positions");
    }
    var geom = {
        type: "LineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// lineString = lineString;
/**
 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
 *
 * @name lineStrings
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
 * @example
 * var linestrings = turf.lineStrings([
 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
 * ]);
 *
 * //=linestrings
 */
function lineStrings(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return lineString(coords, properties);
    }), options);
}
// lineStrings = lineStrings;
/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {FeatureCollection} FeatureCollection of Features
 * @example
 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
 *
 * var collection = turf.featureCollection([
 *   locationA,
 *   locationB,
 *   locationC
 * ]);
 *
 * //=collection
 */
function featureCollection(features, options) {
    if (options === void 0) { options = {}; }
    var fc = { type: "FeatureCollection" };
    if (options.id) {
        fc.id = options.id;
    }
    if (options.bbox) {
        fc.bbox = options.bbox;
    }
    fc.features = features;
    return fc;
}
// featureCollection = featureCollection;
/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 */
function multiLineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiLineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// multiLineString = multiLineString;
/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 */
function multiPoint(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPoint",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// multiPoint = multiPoint;
/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
function multiPolygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPolygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
// multiPolygon = multiPolygon;
/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = turf.geometry("Point", [100, 0]);
 * var line = turf.geometry("LineString", [[101, 0], [102, 1]]);
 * var collection = turf.geometryCollection([pt, line]);
 *
 * // => collection
 */
function geometryCollection(geometries, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "GeometryCollection",
        geometries: geometries,
    };
    return feature(geom, properties, options);
}
// geometryCollection = geometryCollection;
/**
 * Round number to precision
 *
 * @param {number} num Number
 * @param {number} [precision=0] Precision
 * @returns {number} rounded number
 * @example
 * turf.round(120.4321)
 * //=120
 *
 * turf.round(120.4321, 2)
 * //=120.43
 */
function round(num, precision) {
    if (precision === void 0) { precision = 0; }
    if (precision && !(precision >= 0)) {
        throw new Error("precision must be a positive number");
    }
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(num * multiplier) / multiplier;
}
// round = round;
/**
 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name radiansToLength
 * @param {number} radians in radians across the sphere
 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} distance
 */
function radiansToLength(radians, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor =  factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return radians * factor;
}
// radiansToLength = radiansToLength;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name lengthToRadians
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} radians
 */
function lengthToRadians(distance, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor =  factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return distance / factor;
}
// lengthToRadians = lengthToRadians;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
 *
 * @name lengthToDegrees
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} degrees
 */
function lengthToDegrees(distance, units) {
    return radiansToDegrees(lengthToRadians(distance, units));
}
// lengthToDegrees = lengthToDegrees;
/**
 * Converts any bearing angle from the north line direction (positive clockwise)
 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
 *
 * @name bearingToAzimuth
 * @param {number} bearing angle, between -180 and +180 degrees
 * @returns {number} angle between 0 and 360 degrees
 */
function bearingToAzimuth(bearing) {
    var angle = bearing % 360;
    if (angle < 0) {
        angle += 360;
    }
    return angle;
}
// bearingToAzimuth = bearingToAzimuth;
/**
 * Converts an angle in radians to degrees
 *
 * @name radiansToDegrees
 * @param {number} radians angle in radians
 * @returns {number} degrees between 0 and 360 degrees
 */
function radiansToDegrees(radians) {
    var degrees = radians % (2 * Math.PI);
    return (degrees * 180) / Math.PI;
}
// radiansToDegrees = radiansToDegrees;
/**
 * Converts an angle in degrees to radians
 *
 * @name degreesToRadians
 * @param {number} degrees angle between 0 and 360 degrees
 * @returns {number} angle in radians
 */
function degreesToRadians(degrees) {
    var radians = degrees % 360;
    return (radians * Math.PI) / 180;
}
// degreesToRadians = degreesToRadians;
/**
 * Converts a length to the requested unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @param {number} length to be converted
 * @param {Units} [originalUnit="kilometers"] of the length
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted length
 */
function convertLength(length, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "kilometers"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(length >= 0)) {
        throw new Error("length must be a positive number");
    }
    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
}
// convertLength = convertLength;
/**
 * Converts a area to the requested unit.
 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches, hectares
 * @param {number} area to be converted
 * @param {Units} [originalUnit="meters"] of the distance
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted area
 */
function convertArea(area, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "meters"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(area >= 0)) {
        throw new Error("area must be a positive number");
    }
    var startFactor =  areaFactors[originalUnit];
    if (!startFactor) {
        throw new Error("invalid original units");
    }
    var finalFactor =  areaFactors[finalUnit];
    if (!finalFactor) {
        throw new Error("invalid final units");
    }
    return (area / startFactor) * finalFactor;
}
// convertArea = convertArea;
/**
 * isNumber
 *
 * @param {*} num Number to validate
 * @returns {boolean} true/false
 * @example
 * turf.isNumber(123)
 * //=true
 * turf.isNumber('foo')
 * //=false
 */
function isNumber(num) {
    return !isNaN(num) && num !== null && !Array.isArray(num);
}
// isNumber = isNumber;
/**
 * isObject
 *
 * @param {*} input variable to validate
 * @returns {boolean} true/false
 * @example
 * turf.isObject({elevation: 10})
 * //=true
 * turf.isObject('foo')
 * //=false
 */
function isObject(input) {
    return !!input && input.constructor === Object;
}
// isObject = isObject;
/**
 * Validate BBox
 *
 * @private
 * @param {Array<number>} bbox BBox to validate
 * @returns {void}
 * @throws Error if BBox is not valid
 * @example
 * validateBBox([-180, -40, 110, 50])
 * //=OK
 * validateBBox([-180, -40])
 * //=Error
 * validateBBox('Foo')
 * //=Error
 * validateBBox(5)
 * //=Error
 * validateBBox(null)
 * //=Error
 * validateBBox(undefined)
 * //=Error
 */
function validateBBox(bbox) {
    if (!bbox) {
        throw new Error("bbox is required");
    }
    if (!Array.isArray(bbox)) {
        throw new Error("bbox must be an Array");
    }
    if (bbox.length !== 4 && bbox.length !== 6) {
        throw new Error("bbox must be an Array of 4 or 6 numbers");
    }
    bbox.forEach(function (num) {
        if (!isNumber(num)) {
            throw new Error("bbox must only contain numbers");
        }
    });
}
// validateBBox = validateBBox;
/**
 * Validate Id
 *
 * @private
 * @param {string|number} id Id to validate
 * @returns {void}
 * @throws Error if Id is not valid
 * @example
 * validateId([-180, -40, 110, 50])
 * //=Error
 * validateId([-180, -40])
 * //=Error
 * validateId('Foo')
 * //=OK
 * validateId(5)
 * //=OK
 * validateId(null)
 * //=Error
 * validateId(undefined)
 * //=Error
 */
function validateId(id) {
    if (!id) {
        throw new Error("id is required");
    }
    if (["string", "number"].indexOf(typeof id) === -1) {
        throw new Error("id must be a number or a string");
    }
}
// validateId = validateId;















































//'use strict';

//Object.defineProperty(exports, '__esModule', { value: true });

//var helpers = require('@turf/helpers');

/**
 * Callback for coordEach
 *
 * @callback coordEachCallback
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
 *
 * @name coordEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function coordEach(geojson, callback, excludeWrapCoord) {
  // Handles null Geometry -- Skips this GeoJSON
  if (geojson === null) return;
  var j,
    k,
    l,
    geometry,
    stopG,
    coords,
    geometryMaybeCollection,
    wrapShrink = 0,
    coordIndex = 0,
    isGeometryCollection,
    type = geojson.type,
    isFeatureCollection = type === "FeatureCollection",
    isFeature = type === "Feature",
    stop = isFeatureCollection ? geojson.features.length : 1;

  // This logic may look a little weird. The reason why it is that way
  // is because it's trying to be fast. GeoJSON supports multiple kinds
  // of objects at its root: FeatureCollection, Features, Geometries.
  // This function has the responsibility of handling all of them, and that
  // means that some of the `for` loops you see below actually just don't apply
  // to certain inputs. For instance, if you give this just a
  // Point geometry, then both loops are short-circuited and all we do
  // is gradually rename the input until it's called 'geometry'.
  //
  // This also aims to allocate as few resources as possible: just a
  // few numbers and booleans, rather than any temporary arrays as would
  // be required with the normalization approach.
  for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
    geometryMaybeCollection = isFeatureCollection
      ? geojson.features[featureIndex].geometry
      : isFeature
      ? geojson.geometry
      : geojson;
    isGeometryCollection = geometryMaybeCollection
      ? geometryMaybeCollection.type === "GeometryCollection"
      : false;
    stopG = isGeometryCollection
      ? geometryMaybeCollection.geometries.length
      : 1;

    for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
      var multiFeatureIndex = 0;
      var geometryIndex = 0;
      geometry = isGeometryCollection
        ? geometryMaybeCollection.geometries[geomIndex]
        : geometryMaybeCollection;

      // Handles null Geometry -- Skips this geometry
      if (geometry === null) continue;
      coords = geometry.coordinates;
      var geomType = geometry.type;

      wrapShrink =
        excludeWrapCoord &&
        (geomType === "Polygon" || geomType === "MultiPolygon")
          ? 1
          : 0;

      switch (geomType) {
        case null:
          break;
        case "Point":
          if (
            callback(
              coords,
              coordIndex,
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            ) === false
          )
            return false;
          coordIndex++;
          multiFeatureIndex++;
          break;
        case "LineString":
        case "MultiPoint":
          for (j = 0; j < coords.length; j++) {
            if (
              callback(
                coords[j],
                coordIndex,
                featureIndex,
                multiFeatureIndex,
                geometryIndex
              ) === false
            )
              return false;
            coordIndex++;
            if (geomType === "MultiPoint") multiFeatureIndex++;
          }
          if (geomType === "LineString") multiFeatureIndex++;
          break;
        case "Polygon":
        case "MultiLineString":
          for (j = 0; j < coords.length; j++) {
            for (k = 0; k < coords[j].length - wrapShrink; k++) {
              if (
                callback(
                  coords[j][k],
                  coordIndex,
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex
                ) === false
              )
                return false;
              coordIndex++;
            }
            if (geomType === "MultiLineString") multiFeatureIndex++;
            if (geomType === "Polygon") geometryIndex++;
          }
          if (geomType === "Polygon") multiFeatureIndex++;
          break;
        case "MultiPolygon":
          for (j = 0; j < coords.length; j++) {
            geometryIndex = 0;
            for (k = 0; k < coords[j].length; k++) {
              for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                if (
                  callback(
                    coords[j][k][l],
                    coordIndex,
                    featureIndex,
                    multiFeatureIndex,
                    geometryIndex
                  ) === false
                )
                  return false;
                coordIndex++;
              }
              geometryIndex++;
            }
            multiFeatureIndex++;
          }
          break;
        case "GeometryCollection":
          for (j = 0; j < geometry.geometries.length; j++)
            if (
              coordEach(geometry.geometries[j], callback, excludeWrapCoord) ===
              false
            )
              return false;
          break;
        default:
          throw new Error("Unknown Geometry Type");
      }
    }
  }
}






/**
 * Callback for coordReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback coordReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Reduce coordinates in any GeoJSON object, similar to Array.reduce()
 *
 * @name coordReduce
 * @param {FeatureCollection|Geometry|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentCoord, coordIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordReduce(features, function (previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentCoord;
 * });
 */
function coordReduce(geojson, callback, initialValue, excludeWrapCoord) {
  var previousValue = initialValue;
  coordEach(
    geojson,
    function (
      currentCoord,
      coordIndex,
      featureIndex,
      multiFeatureIndex,
      geometryIndex
    ) {
      if (coordIndex === 0 && initialValue === undefined)
        previousValue = currentCoord;
      else
        previousValue = callback(
          previousValue,
          currentCoord,
          coordIndex,
          featureIndex,
          multiFeatureIndex,
          geometryIndex
        );
    },
    excludeWrapCoord
  );
  return previousValue;
}

/**
 * Callback for propEach
 *
 * @callback propEachCallback
 * @param {Object} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over properties in any GeoJSON object, similar to Array.forEach()
 *
 * @name propEach
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentProperties, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propEach(features, function (currentProperties, featureIndex) {
 *   //=currentProperties
 *   //=featureIndex
 * });
 */
function propEach(geojson, callback) {
  var i;
  switch (geojson.type) {
    case "FeatureCollection":
      for (i = 0; i < geojson.features.length; i++) {
        if (callback(geojson.features[i].properties, i) === false) break;
      }
      break;
    case "Feature":
      callback(geojson.properties, 0);
      break;
  }
}

/**
 * Callback for propReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback propReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {*} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce properties in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all properties is unnecessary.
 *
 * @name propReduce
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentProperties, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propReduce(features, function (previousValue, currentProperties, featureIndex) {
 *   //=previousValue
 *   //=currentProperties
 *   //=featureIndex
 *   return currentProperties
 * });
 */
function propReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  propEach(geojson, function (currentProperties, featureIndex) {
    if (featureIndex === 0 && initialValue === undefined)
      previousValue = currentProperties;
    else
      previousValue = callback(previousValue, currentProperties, featureIndex);
  });
  return previousValue;
}

/**
 * Callback for featureEach
 *
 * @callback featureEachCallback
 * @param {Feature<any>} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name featureEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.featureEach(features, function (currentFeature, featureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 * });
 */
function featureEach(geojson, callback) {
  if (geojson.type === "Feature") {
    callback(geojson, 0);
  } else if (geojson.type === "FeatureCollection") {
    for (var i = 0; i < geojson.features.length; i++) {
      if (callback(geojson.features[i], i) === false) break;
    }
  }
}

/**
 * Callback for featureReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback featureReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name featureReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.featureReduce(features, function (previousValue, currentFeature, featureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   return currentFeature
 * });
 */
function featureReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  featureEach(geojson, function (currentFeature, featureIndex) {
    if (featureIndex === 0 && initialValue === undefined)
      previousValue = currentFeature;
    else previousValue = callback(previousValue, currentFeature, featureIndex);
  });
  return previousValue;
}

/**
 * Get all coordinates from any GeoJSON object.
 *
 * @name coordAll
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @returns {Array<Array<number>>} coordinate position array
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * var coords = turf.coordAll(features);
 * //= [[26, 37], [36, 53]]
 */
function coordAll(geojson) {
  var coords = [];
  coordEach(geojson, function (coord) {
    coords.push(coord);
  });
  return coords;
}

/**
 * Callback for geomEach
 *
 * @callback geomEachCallback
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
 *
 * @name geomEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomEach(features, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 * });
 */
function geomEach(geojson, callback) {
  var i,
    j,
    g,
    geometry,
    stopG,
    geometryMaybeCollection,
    isGeometryCollection,
    featureProperties,
    featureBBox,
    featureId,
    featureIndex = 0,
    isFeatureCollection = geojson.type === "FeatureCollection",
    isFeature = geojson.type === "Feature",
    stop = isFeatureCollection ? geojson.features.length : 1;

  // This logic may look a little weird. The reason why it is that way
  // is because it's trying to be fast. GeoJSON supports multiple kinds
  // of objects at its root: FeatureCollection, Features, Geometries.
  // This function has the responsibility of handling all of them, and that
  // means that some of the `for` loops you see below actually just don't apply
  // to certain inputs. For instance, if you give this just a
  // Point geometry, then both loops are short-circuited and all we do
  // is gradually rename the input until it's called 'geometry'.
  //
  // This also aims to allocate as few resources as possible: just a
  // few numbers and booleans, rather than any temporary arrays as would
  // be required with the normalization approach.
  for (i = 0; i < stop; i++) {
    geometryMaybeCollection = isFeatureCollection
      ? geojson.features[i].geometry
      : isFeature
      ? geojson.geometry
      : geojson;
    featureProperties = isFeatureCollection
      ? geojson.features[i].properties
      : isFeature
      ? geojson.properties
      : {};
    featureBBox = isFeatureCollection
      ? geojson.features[i].bbox
      : isFeature
      ? geojson.bbox
      : undefined;
    featureId = isFeatureCollection
      ? geojson.features[i].id
      : isFeature
      ? geojson.id
      : undefined;
    isGeometryCollection = geometryMaybeCollection
      ? geometryMaybeCollection.type === "GeometryCollection"
      : false;
    stopG = isGeometryCollection
      ? geometryMaybeCollection.geometries.length
      : 1;

    for (g = 0; g < stopG; g++) {
      geometry = isGeometryCollection
        ? geometryMaybeCollection.geometries[g]
        : geometryMaybeCollection;

      // Handle null Geometry
      if (geometry === null) {
        if (
          callback(
            null,
            featureIndex,
            featureProperties,
            featureBBox,
            featureId
          ) === false
        )
          return false;
        continue;
      }
      switch (geometry.type) {
        case "Point":
        case "LineString":
        case "MultiPoint":
        case "Polygon":
        case "MultiLineString":
        case "MultiPolygon": {
          if (
            callback(
              geometry,
              featureIndex,
              featureProperties,
              featureBBox,
              featureId
            ) === false
          )
            return false;
          break;
        }
        case "GeometryCollection": {
          for (j = 0; j < geometry.geometries.length; j++) {
            if (
              callback(
                geometry.geometries[j],
                featureIndex,
                featureProperties,
                featureBBox,
                featureId
              ) === false
            )
              return false;
          }
          break;
        }
        default:
          throw new Error("Unknown Geometry Type");
      }
    }
    // Only increase `featureIndex` per each feature
    featureIndex++;
  }
}

/**
 * Callback for geomReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback geomReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Reduce geometry in any GeoJSON object, similar to Array.reduce().
 *
 * @name geomReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomReduce(features, function (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=previousValue
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 *   return currentGeometry
 * });
 */
function geomReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  geomEach(
    geojson,
    function (
      currentGeometry,
      featureIndex,
      featureProperties,
      featureBBox,
      featureId
    ) {
      if (featureIndex === 0 && initialValue === undefined)
        previousValue = currentGeometry;
      else
        previousValue = callback(
          previousValue,
          currentGeometry,
          featureIndex,
          featureProperties,
          featureBBox,
          featureId
        );
    }
  );
  return previousValue;
}

/**
 * Callback for flattenEach
 *
 * @callback flattenEachCallback
 * @param {Feature} currentFeature The current flattened feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Iterate over flattened features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name flattenEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex, multiFeatureIndex)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenEach(features, function (currentFeature, featureIndex, multiFeatureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 * });
 */
function flattenEach(geojson, callback) {
  geomEach(geojson, function (geometry, featureIndex, properties, bbox, id) {
    // Callback for single geometry
    var type = geometry === null ? null : geometry.type;
    switch (type) {
      case null:
      case "Point":
      case "LineString":
      case "Polygon":
        if (
          callback(
            helperfeature(geometry, properties, { bbox: bbox, id: id }),
            featureIndex,
            0
          ) === false
        )
          return false;
        return;
    }

    var geomType;

    // Callback for multi-geometry
    switch (type) {
      case "MultiPoint":
        geomType = "Point";
        break;
      case "MultiLineString":
        geomType = "LineString";
        break;
      case "MultiPolygon":
        geomType = "Polygon";
        break;
    }

    for (
      var multiFeatureIndex = 0;
      multiFeatureIndex < geometry.coordinates.length;
      multiFeatureIndex++
    ) {
      var coordinate = geometry.coordinates[multiFeatureIndex];
      var geom = {
        type: geomType,
        coordinates: coordinate,
      };
      if (
        callback(helperfeature(geom, properties), featureIndex, multiFeatureIndex) ===
        false
      )
        return false;
    }
  });
}

/**
 * Callback for flattenReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback flattenReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Reduce flattened features in any GeoJSON object, similar to Array.reduce().
 *
 * @name flattenReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex, multiFeatureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenReduce(features, function (previousValue, currentFeature, featureIndex, multiFeatureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   return currentFeature
 * });
 */
function flattenReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  flattenEach(
    geojson,
    function (currentFeature, featureIndex, multiFeatureIndex) {
      if (
        featureIndex === 0 &&
        multiFeatureIndex === 0 &&
        initialValue === undefined
      )
        previousValue = currentFeature;
      else
        previousValue = callback(
          previousValue,
          currentFeature,
          featureIndex,
          multiFeatureIndex
        );
    }
  );
  return previousValue;
}

/**
 * Callback for segmentEach
 *
 * @callback segmentEachCallback
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 * @returns {void}
 */

/**
 * Iterate over 2-vertex line segment in any GeoJSON object, similar to Array.forEach()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex)
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentEach(polygon, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //=currentSegment
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   //=segmentIndex
 * });
 *
 * // Calculate the total number of segments
 * var total = 0;
 * turf.segmentEach(polygon, function () {
 *     total++;
 * });
 */
function segmentEach(geojson, callback) {
  flattenEach(geojson, function (feature, featureIndex, multiFeatureIndex) {
    var segmentIndex = 0;

    // Exclude null Geometries
    if (!feature.geometry) return;
    // (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
    var type = feature.geometry.type;
    if (type === "Point" || type === "MultiPoint") return;

    // Generate 2-vertex line segments
    var previousCoords;
    var previousFeatureIndex = 0;
    var previousMultiIndex = 0;
    var prevGeomIndex = 0;
    if (
      coordEach(
        feature,
        function (
          currentCoord,
          coordIndex,
          featureIndexCoord,
          multiPartIndexCoord,
          geometryIndex
        ) {
          // Simulating a meta.coordReduce() since `reduce` operations cannot be stopped by returning `false`
          if (
            previousCoords === undefined ||
            featureIndex > previousFeatureIndex ||
            multiPartIndexCoord > previousMultiIndex ||
            geometryIndex > prevGeomIndex
          ) {
            previousCoords = currentCoord;
            previousFeatureIndex = featureIndex;
            previousMultiIndex = multiPartIndexCoord;
            prevGeomIndex = geometryIndex;
            segmentIndex = 0;
            return;
          }
          var currentSegment = lineString(
            [previousCoords, currentCoord],
            feature.properties
          );
          if (
            callback(
              currentSegment,
              featureIndex,
              multiFeatureIndex,
              geometryIndex,
              segmentIndex
            ) === false
          )
            return false;
          segmentIndex++;
          previousCoords = currentCoord;
        }
      ) === false
    )
      return false;
  });
}

/**
 * Callback for segmentReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback segmentReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 */

/**
 * Reduce 2-vertex line segment in any GeoJSON object, similar to Array.reduce()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (previousValue, currentSegment, currentIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentReduce(polygon, function (previousSegment, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //= previousSegment
 *   //= currentSegment
 *   //= featureIndex
 *   //= multiFeatureIndex
 *   //= geometryIndex
 *   //= segmentIndex
 *   return currentSegment
 * });
 *
 * // Calculate the total number of segments
 * var initialValue = 0
 * var total = turf.segmentReduce(polygon, function (previousValue) {
 *     previousValue++;
 *     return previousValue;
 * }, initialValue);
 */
function segmentReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  var started = false;
  segmentEach(
    geojson,
    function (
      currentSegment,
      featureIndex,
      multiFeatureIndex,
      geometryIndex,
      segmentIndex
    ) {
      if (started === false && initialValue === undefined)
        previousValue = currentSegment;
      else
        previousValue = callback(
          previousValue,
          currentSegment,
          featureIndex,
          multiFeatureIndex,
          geometryIndex,
          segmentIndex
        );
      started = true;
    }
  );
  return previousValue;
}

/**
 * Callback for lineEach
 *
 * @callback lineEachCallback
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Iterate over line or ring coordinates in LineString, Polygon, MultiLineString, MultiPolygon Features or Geometries,
 * similar to Array.forEach.
 *
 * @name lineEach
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @example
 * var multiLine = turf.multiLineString([
 *   [[26, 37], [35, 45]],
 *   [[36, 53], [38, 50], [41, 55]]
 * ]);
 *
 * turf.lineEach(multiLine, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function lineEach(geojson, callback) {
  // validation
  if (!geojson) throw new Error("geojson is required");

  flattenEach(geojson, function (feature, featureIndex, multiFeatureIndex) {
    if (feature.geometry === null) return;
    var type = feature.geometry.type;
    var coords = feature.geometry.coordinates;
    switch (type) {
      case "LineString":
        if (callback(feature, featureIndex, multiFeatureIndex, 0, 0) === false)
          return false;
        break;
      case "Polygon":
        for (
          var geometryIndex = 0;
          geometryIndex < coords.length;
          geometryIndex++
        ) {
          if (
            callback(
              helpers.lineString(coords[geometryIndex], feature.properties),
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            ) === false
          )
            return false;
        }
        break;
    }
  });
}

/**
 * Callback for lineReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback lineReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name lineReduce
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var multiPoly = turf.multiPolygon([
 *   turf.polygon([[[12,48],[2,41],[24,38],[12,48]], [[9,44],[13,41],[13,45],[9,44]]]),
 *   turf.polygon([[[5, 5], [0, 0], [2, 2], [4, 4], [5, 5]]])
 * ]);
 *
 * turf.lineReduce(multiPoly, function (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentLine
 * });
 */
function lineReduce(geojson, callback, initialValue) {
  var previousValue = initialValue;
  lineEach(
    geojson,
    function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
      if (featureIndex === 0 && initialValue === undefined)
        previousValue = currentLine;
      else
        previousValue = callback(
          previousValue,
          currentLine,
          featureIndex,
          multiFeatureIndex,
          geometryIndex
        );
    }
  );
  return previousValue;
}

/**
 * Finds a particular 2-vertex LineString Segment from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 * Point & MultiPoint will always return null.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.segmentIndex=0] Segment Index
 * @param {Object} [options.properties={}] Translate Properties to output LineString
 * @param {BBox} [options.bbox={}] Translate BBox to output LineString
 * @param {number|string} [options.id={}] Translate Id to output LineString
 * @returns {Feature<LineString>} 2-vertex GeoJSON Feature LineString
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findSegment(multiLine);
 * // => Feature<LineString<[[10, 10], [50, 30]]>>
 *
 * // First Segment of 2nd Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: 1});
 * // => Feature<LineString<[[-10, -10], [-50, -30]]>>
 *
 * // Last Segment of Last Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: -1, segmentIndex: -1});
 * // => Feature<LineString<[[-50, -30], [-30, -40]]>>
 */
function findSegment(geojson, options) {
  // Optional Parameters
  options = options || {};
  if (!helpers.isObject(options)) throw new Error("options is invalid");
  var featureIndex = options.featureIndex || 0;
  var multiFeatureIndex = options.multiFeatureIndex || 0;
  var geometryIndex = options.geometryIndex || 0;
  var segmentIndex = options.segmentIndex || 0;

  // Find FeatureIndex
  var properties = options.properties;
  var geometry;

  switch (geojson.type) {
    case "FeatureCollection":
      if (featureIndex < 0)
        featureIndex = geojson.features.length + featureIndex;
      properties = properties || geojson.features[featureIndex].properties;
      geometry = geojson.features[featureIndex].geometry;
      break;
    case "Feature":
      properties = properties || geojson.properties;
      geometry = geojson.geometry;
      break;
    case "Point":
    case "MultiPoint":
      return null;
    case "LineString":
    case "Polygon":
    case "MultiLineString":
    case "MultiPolygon":
      geometry = geojson;
      break;
    default:
      throw new Error("geojson is invalid");
  }

  // Find SegmentIndex
  if (geometry === null) return null;
  var coords = geometry.coordinates;
  switch (geometry.type) {
    case "Point":
    case "MultiPoint":
      return null;
    case "LineString":
      if (segmentIndex < 0) segmentIndex = coords.length + segmentIndex - 1;
      return helpers.lineString(
        [coords[segmentIndex], coords[segmentIndex + 1]],
        properties,
        options
      );
    case "Polygon":
      if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
      if (segmentIndex < 0)
        segmentIndex = coords[geometryIndex].length + segmentIndex - 1;
      return helpers.lineString(
        [
          coords[geometryIndex][segmentIndex],
          coords[geometryIndex][segmentIndex + 1],
        ],
        properties,
        options
      );
    case "MultiLineString":
      if (multiFeatureIndex < 0)
        multiFeatureIndex = coords.length + multiFeatureIndex;
      if (segmentIndex < 0)
        segmentIndex = coords[multiFeatureIndex].length + segmentIndex - 1;
      return helpers.lineString(
        [
          coords[multiFeatureIndex][segmentIndex],
          coords[multiFeatureIndex][segmentIndex + 1],
        ],
        properties,
        options
      );
    case "MultiPolygon":
      if (multiFeatureIndex < 0)
        multiFeatureIndex = coords.length + multiFeatureIndex;
      if (geometryIndex < 0)
        geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
      if (segmentIndex < 0)
        segmentIndex =
          coords[multiFeatureIndex][geometryIndex].length - segmentIndex - 1;
      return helpers.lineString(
        [
          coords[multiFeatureIndex][geometryIndex][segmentIndex],
          coords[multiFeatureIndex][geometryIndex][segmentIndex + 1],
        ],
        properties,
        options
      );
  }
  throw new Error("geojson is invalid");
}

/**
 * Finds a particular Point from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.coordIndex=0] Coord Index
 * @param {Object} [options.properties={}] Translate Properties to output Point
 * @param {BBox} [options.bbox={}] Translate BBox to output Point
 * @param {number|string} [options.id={}] Translate Id to output Point
 * @returns {Feature<Point>} 2-vertex GeoJSON Feature Point
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findPoint(multiLine);
 * // => Feature<Point<[10, 10]>>
 *
 * // First Segment of the 2nd Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: 1});
 * // => Feature<Point<[-10, -10]>>
 *
 * // Last Segment of last Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: -1, coordIndex: -1});
 * // => Feature<Point<[-30, -40]>>
 */
function findPoint(geojson, options) {
  // Optional Parameters
  options = options || {};
  if (!helpers.isObject(options)) throw new Error("options is invalid");
  var featureIndex = options.featureIndex || 0;
  var multiFeatureIndex = options.multiFeatureIndex || 0;
  var geometryIndex = options.geometryIndex || 0;
  var coordIndex = options.coordIndex || 0;

  // Find FeatureIndex
  var properties = options.properties;
  var geometry;

  switch (geojson.type) {
    case "FeatureCollection":
      if (featureIndex < 0)
        featureIndex = geojson.features.length + featureIndex;
      properties = properties || geojson.features[featureIndex].properties;
      geometry = geojson.features[featureIndex].geometry;
      break;
    case "Feature":
      properties = properties || geojson.properties;
      geometry = geojson.geometry;
      break;
    case "Point":
    case "MultiPoint":
      return null;
    case "LineString":
    case "Polygon":
    case "MultiLineString":
    case "MultiPolygon":
      geometry = geojson;
      break;
    default:
      throw new Error("geojson is invalid");
  }

  // Find Coord Index
  if (geometry === null) return null;
  var coords = geometry.coordinates;
  switch (geometry.type) {
    case "Point":
      return helpers.point(coords, properties, options);
    case "MultiPoint":
      if (multiFeatureIndex < 0)
        multiFeatureIndex = coords.length + multiFeatureIndex;
      return helpers.point(coords[multiFeatureIndex], properties, options);
    case "LineString":
      if (coordIndex < 0) coordIndex = coords.length + coordIndex;
      return helpers.point(coords[coordIndex], properties, options);
    case "Polygon":
      if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
      if (coordIndex < 0)
        coordIndex = coords[geometryIndex].length + coordIndex;
      return helpers.point(coords[geometryIndex][coordIndex], properties, options);
    case "MultiLineString":
      if (multiFeatureIndex < 0)
        multiFeatureIndex = coords.length + multiFeatureIndex;
      if (coordIndex < 0)
        coordIndex = coords[multiFeatureIndex].length + coordIndex;
      return helpers.point(coords[multiFeatureIndex][coordIndex], properties, options);
    case "MultiPolygon":
      if (multiFeatureIndex < 0)
        multiFeatureIndex = coords.length + multiFeatureIndex;
      if (geometryIndex < 0)
        geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
      if (coordIndex < 0)
        coordIndex =
          coords[multiFeatureIndex][geometryIndex].length - coordIndex;
      return helpers.point(
        coords[multiFeatureIndex][geometryIndex][coordIndex],
        properties,
        options
      );
  }
  throw new Error("geojson is invalid");
}

 coordAll = coordAll;
 coordEach = coordEach;
 coordReduce = coordReduce;
 featureEach = featureEach;
 featureReduce = featureReduce;
 findPoint = findPoint;
 findSegment = findSegment;
 flattenEach = flattenEach;
 flattenReduce = flattenReduce;
 geomEach = geomEach;
 geomReduce = geomReduce;
 lineEach = lineEach;
 lineReduce = lineReduce;
 propEach = propEach;
 propReduce = propReduce;
 segmentEach = segmentEach;
 segmentReduce = segmentReduce;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 //"use strict";
//Object.defineProperty(exports, "__esModule", { value: true });
//var invariant_1 = require("@turf/invariant");
//var helpers_1 = require("@turf/helpers");
//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html
/**
 * Calculates the distance between two {@link Point|points} in degrees, radians, miles, or kilometers.
 * This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
 *
 * @name distance
 * @param {Coord | Point} from origin point or coordinate
 * @param {Coord | Point} to destination point or coordinate
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
 * @returns {number} distance between the two points
 * @example
 * var from = turf.point([-75.343, 39.984]);
 * var to = turf.point([-75.534, 39.123]);
 * var options = {units: 'miles'};
 *
 * var distance = turf.distance(from, to, options);
 *
 * //addToMap
 * var addToMap = [from, to];
 * from.properties.distance = distance;
 * to.properties.distance = distance;
 */
function distance(from, to, options) {
    if (options === void 0) { options = {}; }
    var coordinates1 = getCoord(from);
    var coordinates2 = getCoord(to);
    var dLat = degreesToRadians(coordinates2[1] - coordinates1[1]);
    var dLon = degreesToRadians(coordinates2[0] - coordinates1[0]);
    var lat1 = degreesToRadians(coordinates1[1]);
    var lat2 = degreesToRadians(coordinates2[1]);
    var a = Math.pow(Math.sin(dLat / 2), 2) +
        Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
    return radiansToLength(2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)), options.units);
}

//exports.default = distance;



















////"use strict";
//Object.defineProperty(exports, "__esModule", { value: true });
//var helpers_1 = require("@turf/helpers");
/**
 * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
 *
 * @name getCoord
 * @param {Array<number>|Geometry<Point>|Feature<Point>} coord GeoJSON Point or an Array of numbers
 * @returns {Array<number>} coordinates
 * @example
 * var pt = turf.point([10, 10]);
 *
 * var coord = turf.getCoord(pt);
 * //= [10, 10]
 */
 
 
 
function getCoord(coord) {
    if (!coord) {
        throw new Error("coord is required");
    }
    if (!Array.isArray(coord)) {
        if (coord.type === "Feature" &&
            coord.geometry !== null &&
            coord.geometry.type === "Point") {
            return coord.geometry.coordinates;
        }
        if (coord.type === "Point") {
            return coord.coordinates;
        }
    }
    if (Array.isArray(coord) &&
        coord.length >= 2 &&
        !Array.isArray(coord[0]) &&
        !Array.isArray(coord[1])) {
        return coord;
    }
    throw new Error("coord must be GeoJSON Point or an Array of numbers");
}
//exports.getCoord = getCoord;
/**
 * Unwrap coordinates from a Feature, Geometry Object or an Array
 *
 * @name getCoords
 * @param {Array<any>|Geometry|Feature} coords Feature, Geometry Object or an Array
 * @returns {Array<any>} coordinates
 * @example
 * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
 *
 * var coords = turf.getCoords(poly);
 * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
 */
function getCoords(coords) {
    if (Array.isArray(coords)) {
        return coords;
    }
    // Feature
    if (coords.type === "Feature") {
        if (coords.geometry !== null) {
            return coords.geometry.coordinates;
        }
    }
    else {
        // Geometry
        if (coords.coordinates) {
            return coords.coordinates;
        }
    }
    throw new Error("coords must be GeoJSON Feature, Geometry Object or an Array");
}
//exports.getCoords = getCoords;
/**
 * Checks if coordinates contains a number
 *
 * @name containsNumber
 * @param {Array<any>} coordinates GeoJSON Coordinates
 * @returns {boolean} true if Array contains a number
 */
function containsNumber(coordinates) {
    if (coordinates.length > 1 &&
        isNumber(coordinates[0]) &&
        isNumber(coordinates[1])) {
        return true;
    }
    if (Array.isArray(coordinates[0]) && coordinates[0].length) {
        return containsNumber(coordinates[0]);
    }
    throw new Error("coordinates must only contain numbers");
}
//exports.containsNumber = containsNumber;
/**
 * Enforce expectations about types of GeoJSON objects for Turf.
 *
 * @name geojsonType
 * @param {GeoJSON} value any GeoJSON object
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function geojsonType(value, type, name) {
    if (!type || !name) {
        throw new Error("type and name required");
    }
    if (!value || value.type !== type) {
        throw new Error("Invalid input to " +
            name +
            ": must be a " +
            type +
            ", given " +
            value.type);
    }
}
//exports.geojsonType = geojsonType;
/**
 * Enforce expectations about types of {@link Feature} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name featureOf
 * @param {Feature} feature a feature with an expected geometry type
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} error if value is not the expected type.
 */
function featureOf(feature, type, name) {
    if (!feature) {
        throw new Error("No feature passed");
    }
    if (!name) {
        throw new Error(".featureOf() requires a name");
    }
    if (!feature || feature.type !== "Feature" || !feature.geometry) {
        throw new Error("Invalid input to " + name + ", Feature with geometry required");
    }
    if (!feature.geometry || feature.geometry.type !== type) {
        throw new Error("Invalid input to " +
            name +
            ": must be a " +
            type +
            ", given " +
            feature.geometry.type);
    }
}
//exports.featureOf = featureOf;
/**
 * Enforce expectations about types of {@link FeatureCollection} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name collectionOf
 * @param {FeatureCollection} featureCollection a FeatureCollection for which features will be judged
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function collectionOf(featureCollection, type, name) {
    if (!featureCollection) {
        throw new Error("No featureCollection passed");
    }
    if (!name) {
        throw new Error(".collectionOf() requires a name");
    }
    if (!featureCollection || featureCollection.type !== "FeatureCollection") {
        throw new Error("Invalid input to " + name + ", FeatureCollection required");
    }
    for (var _i = 0, _a = featureCollection.features; _i < _a.length; _i++) {
        var feature = _a[_i];
        if (!feature || feature.type !== "Feature" || !feature.geometry) {
            throw new Error("Invalid input to " + name + ", Feature with geometry required");
        }
        if (!feature.geometry || feature.geometry.type !== type) {
            throw new Error("Invalid input to " +
                name +
                ": must be a " +
                type +
                ", given " +
                feature.geometry.type);
        }
    }
}
//exports.collectionOf = collectionOf;
/**
 * Get Geometry from Feature or Geometry Object
 *
 * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
 * @returns {Geometry|null} GeoJSON Geometry Object
 * @throws {Error} if geojson is not a Feature or Geometry Object
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getGeom(point)
 * //={"type": "Point", "coordinates": [110, 40]}
 */
function getGeom(geojson) {
    if (geojson.type === "Feature") {
        return geojson.geometry;
    }
    return geojson;
}
//exports.getGeom = getGeom;
/**
 * Get GeoJSON object's type, Geometry type is prioritize.
 *
 * @param {GeoJSON} geojson GeoJSON object
 * @param {string} [name="geojson"] name of the variable to display in error message (unused)
 * @returns {string} GeoJSON type
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getType(point)
 * //="Point"
 */
function getType(geojson, _name) {
    if (geojson.type === "FeatureCollection") {
        return "FeatureCollection";
    }
    if (geojson.type === "GeometryCollection") {
        return "GeometryCollection";
    }
    if (geojson.type === "Feature" && geojson.geometry !== null) {
        return geojson.geometry.type;
    }
    return geojson.type;
}
//exports.getType = getType;























var each = coordEach;

/**
 * Takes any {@link GeoJSON} object, calculates the extent of all input features, and returns a bounding box.
 *
 * @module turf/extent
 * @category measurement
 * @param {GeoJSON} input any valid GeoJSON Object
 * @return {Array<number>} the bounding box of `input` given
 * as an array in WSEN order (west, south, east, north)
 * @example
 * var input = {
 *   "type": "FeatureCollection",
 *   "features": [
 *     {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [114.175329, 22.2524]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [114.170007, 22.267969]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [114.200649, 22.274641]
 *       }
 *     }, {
 *       "type": "Feature",
 *       "properties": {},
 *       "geometry": {
 *         "type": "Point",
 *         "coordinates": [114.186744, 22.265745]
 *       }
 *     }
 *   ]
 * };
 *
 * var bbox = turf.extent(input);
 *
 * var bboxPolygon = turf.bboxPolygon(bbox);
 *
 * var resultFeatures = input.features.concat(bboxPolygon);
 * var result = {
 *   "type": "FeatureCollection",
 *   "features": resultFeatures
 * };
 *
 * //=result
 */
//module.exports 
 function extent(layer) {
    var extent = [Infinity, Infinity, -Infinity, -Infinity];
    each(layer, function(coord) {
      if (extent[0] > coord[0]) extent[0] = coord[0];
      if (extent[1] > coord[1]) extent[1] = coord[1];
      if (extent[2] < coord[0]) extent[2] = coord[0];
      if (extent[3] < coord[1]) extent[3] = coord[1];
    });
    return extent;
}





























































































































//start routing

(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(_dereq_,module,exports){
function corslite(url, callback, cors) {
    var sent = false;

    if (typeof window.XMLHttpRequest === 'undefined') {
        return callback(Error('Browser not supported'));
    }

    if (typeof cors === 'undefined') {
        var m = url.match(/^\s*https?:\/\/[^\/]*/);
        cors = m && (m[0] !== location.protocol + '//' + location.hostname +
                (location.port ? ':' + location.port : ''));
    }

    var x = new window.XMLHttpRequest();

    function isSuccessful(status) {
        return status >= 200 && status < 300 || status === 304;
    }

    if (cors && !('withCredentials' in x)) {
        // IE8-9
        x = new window.XDomainRequest();

        // Ensure callback is never called synchronously, i.e., before
        // x.send() returns (this has been observed in the wild).
        // See https://github.com/mapbox/mapbox.js/issues/472
        var original = callback;
        callback = function() {
            if (sent) {
                original.apply(this, arguments);
            } else {
                var that = this, args = arguments;
                setTimeout(function() {
                    original.apply(that, args);
                }, 0);
            }
        }
    }

    function loaded() {
        if (
            // XDomainRequest
            x.status === undefined ||
            // modern browsers
            isSuccessful(x.status)) callback.call(x, null, x);
        else callback.call(x, x, null);
    }

    // Both `onreadystatechange` and `onload` can fire. `onreadystatechange`
    // has [been supported for longer](http://stackoverflow.com/a/9181508/229001).
    if ('onload' in x) {
        x.onload = loaded;
    } else {
        x.onreadystatechange = function readystate() {
            if (x.readyState === 4) {
                loaded();
            }
        };
    }

    // Call the callback with the XMLHttpRequest object as an error and prevent
    // it from ever being called again by reassigning it to `noop`
    x.onerror = function error(evt) {
        // XDomainRequest provides no evt parameter
        callback.call(this, evt || true, null);
        callback = function() { };
    };

    // IE9 must have onprogress be set to a unique function.
    x.onprogress = function() { };

    x.ontimeout = function(evt) {
        callback.call(this, evt, null);
        callback = function() { };
    };

    x.onabort = function(evt) {
        callback.call(this, evt, null);
        callback = function() { };
    };

    // GET is the only supported HTTP Verb by XDomainRequest and is the
    // only one supported here.
    x.open('GET', url, true);

    // Send the request. Sending data is not supported.
    x.send(null);
    sent = true;

    return x;
}

if (typeof module !== 'undefined') module.exports = corslite;

},{}],2:[function(_dereq_,module,exports){
//'use strict';

/**
 * Based off of [the offical Google document](https://developers.google.com/maps/documentation/utilities/polylinealgorithm)
 *
 * Some parts from [this implementation](http://facstaff.unca.edu/mcmcclur/GoogleMaps/EncodePolyline/PolylineEncoder.js)
 * by [Mark McClure](http://facstaff.unca.edu/mcmcclur/)
 *
 * @module polyline
 */

var polyline = {};

function py2_round(value) {
    // Google's polyline algorithm uses the same rounding strategy as Python 2, which is different from JS for negative values
    return Math.floor(Math.abs(value) + 0.5) * Math.sign(value);
}

function encode(current, previous, factor) {
    current = py2_round(current * factor);
    previous = py2_round(previous * factor);
    var coordinate = current - previous;
    coordinate <<= 1;
    if (current - previous < 0) {
        coordinate = ~coordinate;
    }
    var output = '';
    while (coordinate >= 0x20) {
        output += String.fromCharCode((0x20 | (coordinate & 0x1f)) + 63);
        coordinate >>= 5;
    }
    output += String.fromCharCode(coordinate + 63);
    return output;
}

/**
 * Decodes to a [latitude, longitude] coordinates array.
 *
 * This is adapted from the implementation in Project-OSRM.
 *
 * @param {String} str
 * @param {Number} precision
 * @returns {Array}
 *
 * @see https://github.com/Project-OSRM/osrm-frontend/blob/master/WebContent/routing/OSRM.RoutingGeometry.js
 */
polyline.decode = function(str, precision) {
    var index = 0,
        lat = 0,
        lng = 0,
        coordinates = [],
        shift = 0,
        result = 0,
        byte = null,
        latitude_change,
        longitude_change,
        factor = Math.pow(10, precision || 5);

    // Coordinates have variable length when encoded, so just keep
    // track of whether we've hit the end of the string. In each
    // loop iteration, a single coordinate is decoded.
    while (index < str.length) {

        // Reset shift, result, and byte
        byte = null;
        shift = 0;
        result = 0;

        do {
            byte = str.charCodeAt(index++) - 63;
            result |= (byte & 0x1f) << shift;
            shift += 5;
        } while (byte >= 0x20);

        latitude_change = ((result & 1) ? ~(result >> 1) : (result >> 1));

        shift = result = 0;

        do {
            byte = str.charCodeAt(index++) - 63;
            result |= (byte & 0x1f) << shift;
            shift += 5;
        } while (byte >= 0x20);

        longitude_change = ((result & 1) ? ~(result >> 1) : (result >> 1));

        lat += latitude_change;
        lng += longitude_change;

        coordinates.push([lat / factor, lng / factor]);
    }

    return coordinates;
};

/**
 * Encodes the given [latitude, longitude] coordinates array.
 *
 * @param {Array.<Array.<Number>>} coordinates
 * @param {Number} precision
 * @returns {String}
 */
polyline.encode = function(coordinates, precision) {
    if (!coordinates.length) { return ''; }

    var factor = Math.pow(10, precision || 5),
        output = encode(coordinates[0][0], 0, factor) + encode(coordinates[0][1], 0, factor);

    for (var i = 1; i < coordinates.length; i++) {
        var a = coordinates[i], b = coordinates[i - 1];
        output += encode(a[0], b[0], factor);
        output += encode(a[1], b[1], factor);
    }

    return output;
};

function flipped(coords) {
    var flipped = [];
    for (var i = 0; i < coords.length; i++) {
        flipped.push(coords[i].slice().reverse());
    }
    return flipped;
}

/**
 * Encodes a GeoJSON LineString feature/geometry.
 *
 * @param {Object} geojson
 * @param {Number} precision
 * @returns {String}
 */
polyline.fromGeoJSON = function(geojson, precision) {
    if (geojson && geojson.type === 'Feature') {
        geojson = geojson.geometry;
    }
    if (!geojson || geojson.type !== 'LineString') {
        throw new Error('Input must be a GeoJSON LineString');
    }
    return polyline.encode(flipped(geojson.coordinates), precision);
};

/**
 * Decodes to a GeoJSON LineString geometry.
 *
 * @param {String} str
 * @param {Number} precision
 * @returns {Object}
 */
polyline.toGeoJSON = function(str, precision) {
    var coords = polyline.decode(str, precision);
    return {
        type: 'LineString',
        coordinates: flipped(coords)
    };
};

if (typeof module === 'object' && module.exports) {
    module.exports = polyline;
}

},{}],3:[function(_dereq_,module,exports){
var languages = _dereq_('./languages');
var instructions = languages.instructions;
var grammars = languages.grammars;
var abbreviations = languages.abbreviations;

module.exports = function(version) {
    Object.keys(instructions).forEach(function(code) {
        if (!instructions[code][version]) { throw 'invalid version ' + version + ': ' + code + ' not supported'; }
    });

    return {
        capitalizeFirstLetter: function(language, string) {
            return string.charAt(0).toLocaleUpperCase(language) + string.slice(1);
        },
        ordinalize: function(language, number) {
            // Transform numbers to their translated ordinalized value
            if (!language) throw new Error('No language code provided');

            return instructions[language][version].constants.ordinalize[number.toString()] || '';
        },
        directionFromDegree: function(language, degree) {
            // Transform degrees to their translated compass direction
            if (!language) throw new Error('No language code provided');
            if (!degree && degree !== 0) {
                // step had no bearing_after degree, ignoring
                return '';
            } else if (degree >= 0 && degree <= 20) {
                return instructions[language][version].constants.direction.north;
            } else if (degree > 20 && degree < 70) {
                return instructions[language][version].constants.direction.northeast;
            } else if (degree >= 70 && degree <= 110) {
                return instructions[language][version].constants.direction.east;
            } else if (degree > 110 && degree < 160) {
                return instructions[language][version].constants.direction.southeast;
            } else if (degree >= 160 && degree <= 200) {
                return instructions[language][version].constants.direction.south;
            } else if (degree > 200 && degree < 250) {
                return instructions[language][version].constants.direction.southwest;
            } else if (degree >= 250 && degree <= 290) {
                return instructions[language][version].constants.direction.west;
            } else if (degree > 290 && degree < 340) {
                return instructions[language][version].constants.direction.northwest;
            } else if (degree >= 340 && degree <= 360) {
                return instructions[language][version].constants.direction.north;
            } else {
                throw new Error('Degree ' + degree + ' invalid');
            }
        },
        laneConfig: function(step) {
            // Reduce any lane combination down to a contracted lane diagram
            if (!step.intersections || !step.intersections[0].lanes) throw new Error('No lanes object');

            var config = [];
            var currentLaneValidity = null;

            step.intersections[0].lanes.forEach(function (lane) {
                if (currentLaneValidity === null || currentLaneValidity !== lane.valid) {
                    if (lane.valid) {
                        config.push('o');
                    } else {
                        config.push('x');
                    }
                    currentLaneValidity = lane.valid;
                }
            });

            return config.join('');
        },
        getWayName: function(language, step, options) {
            var classes = options ? options.classes || [] : [];
            if (typeof step !== 'object') throw new Error('step must be an Object');
            if (!language) throw new Error('No language code provided');
            if (!Array.isArray(classes)) throw new Error('classes must be an Array or undefined');

            var wayName;
            var name = step.name || '';
            var ref = (step.ref || '').split(';')[0];

            // Remove hacks from Mapbox Directions mixing ref into name
            if (name === step.ref) {
                // if both are the same we assume that there used to be an empty name, with the ref being filled in for it
                // we only need to retain the ref then
                name = '';
            }
            name = name.replace(' (' + step.ref + ')', '');

            // In attempt to avoid using the highway name of a way,
            // check and see if the step has a class which should signal
            // the ref should be used instead of the name.
            var wayMotorway = classes.indexOf('motorway') !== -1;

            if (name && ref && name !== ref && !wayMotorway) {
                var phrase = instructions[language][version].phrase['name and ref'] ||
                    instructions.en[version].phrase['name and ref'];
                wayName = this.tokenize(language, phrase, {
                    name: name,
                    ref: ref
                }, options);
            } else if (name && ref && wayMotorway && (/\d/).test(ref)) {
                wayName = options && options.formatToken ? options.formatToken('ref', ref) : ref;
            } else if (!name && ref) {
                wayName = options && options.formatToken ? options.formatToken('ref', ref) : ref;
            } else {
                wayName = options && options.formatToken ? options.formatToken('name', name) : name;
            }

            return wayName;
        },

        /**
         * Formulate a localized text instruction from a step.
         *
         * @param  {string} language           Language code.
         * @param  {object} step               Step including maneuver property.
         * @param  {object} opts               Additional options.
         * @param  {string} opts.legIndex      Index of leg in the route.
         * @param  {string} opts.legCount      Total number of legs in the route.
         * @param  {array}  opts.classes       List of road classes.
         * @param  {string} opts.waypointName  Name of waypoint for arrival instruction.
         *
         * @return {string} Localized text instruction.
         */
        compile: function(language, step, opts) {
            if (!language) throw new Error('No language code provided');
            if (languages.supportedCodes.indexOf(language) === -1) throw new Error('language code ' + language + ' not loaded');
            if (!step.maneuver) throw new Error('No step maneuver provided');
            var options = opts || {};

            var type = step.maneuver.type;
            var modifier = step.maneuver.modifier;
            var mode = step.mode;
            // driving_side will only be defined in OSRM 5.14+
            var side = step.driving_side;

            if (!type) { throw new Error('Missing step maneuver type'); }
            if (type !== 'depart' && type !== 'arrive' && !modifier) { throw new Error('Missing step maneuver modifier'); }

            if (!instructions[language][version][type]) {
                // Log for debugging
                console.log('Encountered unknown instruction type: ' + type); // eslint-disable-line no-console
                // OSRM specification assumes turn types can be added without
                // major version changes. Unknown types are to be treated as
                // type `turn` by clients
                type = 'turn';
            }

            // Use special instructions if available, otherwise `defaultinstruction`
            var instructionObject;
            if (instructions[language][version].modes[mode]) {
                instructionObject = instructions[language][version].modes[mode];
            } else {
              // omit side from off ramp if same as driving_side
              // note: side will be undefined if the input is from OSRM <5.14
              // but the condition should still evaluate properly regardless
                var omitSide = type === 'off ramp' && modifier.indexOf(side) >= 0;
                if (instructions[language][version][type][modifier] && !omitSide) {
                    instructionObject = instructions[language][version][type][modifier];
                } else {
                    instructionObject = instructions[language][version][type].default;
                }
            }

            // Special case handling
            var laneInstruction;
            switch (type) {
            case 'use lane':
                laneInstruction = instructions[language][version].constants.lanes[this.laneConfig(step)];
                if (!laneInstruction) {
                    // If the lane combination is not found, default to continue straight
                    instructionObject = instructions[language][version]['use lane'].no_lanes;
                }
                break;
            case 'rotary':
            case 'roundabout':
                if (step.rotary_name && step.maneuver.exit && instructionObject.name_exit) {
                    instructionObject = instructionObject.name_exit;
                } else if (step.rotary_name && instructionObject.name) {
                    instructionObject = instructionObject.name;
                } else if (step.maneuver.exit && instructionObject.exit) {
                    instructionObject = instructionObject.exit;
                } else {
                    instructionObject = instructionObject.default;
                }
                break;
            default:
                // NOOP, since no special logic for that type
            }

            // Decide way_name with special handling for name and ref
            var wayName = this.getWayName(language, step, options);

            // Decide which instruction string to use
            // Destination takes precedence over name
            var instruction;
            if (step.destinations && step.exits && instructionObject.exit_destination) {
                instruction = instructionObject.exit_destination;
            } else if (step.destinations && instructionObject.destination) {
                instruction = instructionObject.destination;
            } else if (step.exits && instructionObject.exit) {
                instruction = instructionObject.exit;
            } else if (wayName && instructionObject.name) {
                instruction = instructionObject.name;
            } else if (options.waypointName && instructionObject.named) {
                instruction = instructionObject.named;
            } else {
                instruction = instructionObject.default;
            }

            var destinations = step.destinations && step.destinations.split(': ');
            var destinationRef = destinations && destinations[0].split(',')[0];
            var destination = destinations && destinations[1] && destinations[1].split(',')[0];
            var firstDestination;
            if (destination && destinationRef) {
                firstDestination = destinationRef + ': ' + destination;
            } else {
                firstDestination = destinationRef || destination || '';
            }

            var nthWaypoint = options.legIndex >= 0 && options.legIndex !== options.legCount - 1 ? this.ordinalize(language, options.legIndex + 1) : '';

            // Replace tokens
            // NOOP if they don't exist
            var replaceTokens = {
                'way_name': wayName,
                'destination': firstDestination,
                'exit': (step.exits || '').split(';')[0],
                'exit_number': this.ordinalize(language, step.maneuver.exit || 1),
                'rotary_name': step.rotary_name,
                'lane_instruction': laneInstruction,
                'modifier': instructions[language][version].constants.modifier[modifier],
                'direction': this.directionFromDegree(language, step.maneuver.bearing_after),
                'nth': nthWaypoint,
                'waypoint_name': options.waypointName
            };

            return this.tokenize(language, instruction, replaceTokens, options);
        },
        grammarize: function(language, name, grammar) {
            if (!language) throw new Error('No language code provided');
            // Process way/rotary name with applying grammar rules if any
            if (name && grammar && grammars && grammars[language] && grammars[language][version]) {
                var rules = grammars[language][version][grammar];
                if (rules) {
                    // Pass original name to rules' regular expressions enclosed with spaces for simplier parsing
                    var n = ' ' + name + ' ';
                    var flags = grammars[language].meta.regExpFlags || '';
                    rules.forEach(function(rule) {
                        var re = new RegExp(rule[0], flags);
                        n = n.replace(re, rule[1]);
                    });

                    return n.trim();
                }
            }

            return name;
        },
        abbreviations: abbreviations,
        tokenize: function(language, instruction, tokens, options) {
            if (!language) throw new Error('No language code provided');
            // Keep this function context to use in inline function below (no arrow functions in ES4)
            var that = this;
            var startedWithToken = false;
            var output = instruction.replace(/\{(\w+)(?::(\w+))?\}/g, function(token, tag, grammar, offset) {
                var value = tokens[tag];

                // Return unknown token unchanged
                if (typeof value === 'undefined') {
                    return token;
                }

                value = that.grammarize(language, value, grammar);

                // If this token appears at the beginning of the instruction, capitalize it.
                if (offset === 0 && instructions[language].meta.capitalizeFirstLetter) {
                    startedWithToken = true;
                    value = that.capitalizeFirstLetter(language, value);
                }

                if (options && options.formatToken) {
                    value = options.formatToken(tag, value);
                }

                return value;
            })
            .replace(/ {2}/g, ' '); // remove excess spaces

            if (!startedWithToken && instructions[language].meta.capitalizeFirstLetter) {
                return this.capitalizeFirstLetter(language, output);
            }

            return output;
        }
    };
};

},{"./languages":4}],4:[function(_dereq_,module,exports){
// Load all language files explicitly to allow integration
// with bundling tools like webpack and browserify
var instructionsDa = _dereq_('./languages/translations/da.json');
var instructionsDe = _dereq_('./languages/translations/de.json');
var instructionsEn = _dereq_('./languages/translations/en.json');
var instructionsEo = _dereq_('./languages/translations/eo.json');
var instructionsEs = _dereq_('./languages/translations/es.json');
var instructionsEsEs = _dereq_('./languages/translations/es-ES.json');
var instructionsFi = _dereq_('./languages/translations/fi.json');
var instructionsFr = _dereq_('./languages/translations/fr.json');
var instructionsHe = _dereq_('./languages/translations/he.json');
var instructionsId = _dereq_('./languages/translations/id.json');
var instructionsIt = _dereq_('./languages/translations/it.json');
var instructionsKo = _dereq_('./languages/translations/ko.json');
var instructionsMy = _dereq_('./languages/translations/my.json');
var instructionsNl = _dereq_('./languages/translations/nl.json');
var instructionsNo = _dereq_('./languages/translations/no.json');
var instructionsPl = _dereq_('./languages/translations/pl.json');
var instructionsPtBr = _dereq_('./languages/translations/pt-BR.json');
var instructionsPtPt = _dereq_('./languages/translations/pt-PT.json');
var instructionsRo = _dereq_('./languages/translations/ro.json');
var instructionsRu = _dereq_('./languages/translations/ru.json');
var instructionsSv = _dereq_('./languages/translations/sv.json');
var instructionsTr = _dereq_('./languages/translations/tr.json');
var instructionsUk = _dereq_('./languages/translations/uk.json');
var instructionsVi = _dereq_('./languages/translations/vi.json');
var instructionsZhHans = _dereq_('./languages/translations/zh-Hans.json');

// Load all grammar files
var grammarFr = _dereq_('./languages/grammar/fr.json');
var grammarRu = _dereq_('./languages/grammar/ru.json');

// Load all abbreviations files
var abbreviationsBg = _dereq_('./languages/abbreviations/bg.json');
var abbreviationsCa = _dereq_('./languages/abbreviations/ca.json');
var abbreviationsDa = _dereq_('./languages/abbreviations/da.json');
var ebbreviationsDe = _dereq_('./languages/abbreviations/de.json');
var abbreviationsEn = _dereq_('./languages/abbreviations/en.json');
var abbreviationsEs = _dereq_('./languages/abbreviations/es.json');
var abbreviationsFr = _dereq_('./languages/abbreviations/fr.json');
var abbreviationsHe = _dereq_('./languages/abbreviations/he.json');
var abbreviationsHu = _dereq_('./languages/abbreviations/hu.json');
var abbreviationsLt = _dereq_('./languages/abbreviations/lt.json');
var abbreviationsNl = _dereq_('./languages/abbreviations/nl.json');
var abbreviationsRu = _dereq_('./languages/abbreviations/ru.json');
var abbreviationsSl = _dereq_('./languages/abbreviations/sl.json');
var abbreviationsSv = _dereq_('./languages/abbreviations/sv.json');
var abbreviationsUk = _dereq_('./languages/abbreviations/uk.json');
var abbreviationsVi = _dereq_('./languages/abbreviations/vi.json');

// Create a list of supported codes
var instructions = {
    'da': instructionsDa,
    'de': instructionsDe,
    'en': instructionsEn,
    'eo': instructionsEo,
    'es': instructionsEs,
    'es-ES': instructionsEsEs,
    'fi': instructionsFi,
    'fr': instructionsFr,
    'he': instructionsHe,
    'id': instructionsId,
    'it': instructionsIt,
    'ko': instructionsKo,
    'my': instructionsMy,
    'nl': instructionsNl,
    'no': instructionsNo,
    'pl': instructionsPl,
    'pt-BR': instructionsPtBr,
    'pt-PT': instructionsPtPt,
    'ro': instructionsRo,
    'ru': instructionsRu,
    'sv': instructionsSv,
    'tr': instructionsTr,
    'uk': instructionsUk,
    'vi': instructionsVi,
    'zh-Hans': instructionsZhHans
};

// Create list of supported grammar
var grammars = {
    'fr': grammarFr,
    'ru': grammarRu
};

// Create list of supported abbrevations
var abbreviations = {
    'bg': abbreviationsBg,
    'ca': abbreviationsCa,
    'da': abbreviationsDa,
    'de': ebbreviationsDe,
    'en': abbreviationsEn,
    'es': abbreviationsEs,
    'fr': abbreviationsFr,
    'he': abbreviationsHe,
    'hu': abbreviationsHu,
    'lt': abbreviationsLt,
    'nl': abbreviationsNl,
    'ru': abbreviationsRu,
    'sl': abbreviationsSl,
    'sv': abbreviationsSv,
    'uk': abbreviationsUk,
    'vi': abbreviationsVi
};
module.exports = {
    supportedCodes: Object.keys(instructions),
    instructions: instructions,
    grammars: grammars,
    abbreviations: abbreviations
};

},{"./languages/abbreviations/bg.json":5,"./languages/abbreviations/ca.json":6,"./languages/abbreviations/da.json":7,"./languages/abbreviations/de.json":8,"./languages/abbreviations/en.json":9,"./languages/abbreviations/es.json":10,"./languages/abbreviations/fr.json":11,"./languages/abbreviations/he.json":12,"./languages/abbreviations/hu.json":13,"./languages/abbreviations/lt.json":14,"./languages/abbreviations/nl.json":15,"./languages/abbreviations/ru.json":16,"./languages/abbreviations/sl.json":17,"./languages/abbreviations/sv.json":18,"./languages/abbreviations/uk.json":19,"./languages/abbreviations/vi.json":20,"./languages/grammar/fr.json":21,"./languages/grammar/ru.json":22,"./languages/translations/da.json":23,"./languages/translations/de.json":24,"./languages/translations/en.json":25,"./languages/translations/eo.json":26,"./languages/translations/es-ES.json":27,"./languages/translations/es.json":28,"./languages/translations/fi.json":29,"./languages/translations/fr.json":30,"./languages/translations/he.json":31,"./languages/translations/id.json":32,"./languages/translations/it.json":33,"./languages/translations/ko.json":34,"./languages/translations/my.json":35,"./languages/translations/nl.json":36,"./languages/translations/no.json":37,"./languages/translations/pl.json":38,"./languages/translations/pt-BR.json":39,"./languages/translations/pt-PT.json":40,"./languages/translations/ro.json":41,"./languages/translations/ru.json":42,"./languages/translations/sv.json":43,"./languages/translations/tr.json":44,"./languages/translations/uk.json":45,"./languages/translations/vi.json":46,"./languages/translations/zh-Hans.json":47}],5:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "международен": "Межд",
        "старши": "Стрш",
        "възел": "Въз",
        "пазар": "Mkt",
        "светисвети": "СвСв",
        "сестра": "сес",
        "уилям": "Ум",
        "апартаменти": "ап",
        "езеро": "Ез",
        "свети": "Св",
        "център": "Ц-р",
        "парк": "Пк",
        "маршрут": "М-т",
        "площад": "Пл",
        "национален": "Нац",
        "училище": "Уч",
        "река": "Рек",
        "поток": "П-к",
        "район": "Р-н",
        "крепост": "К-т",
        "паметник": "Пам",
        "университет": "Уни",
        "Връх": "Вр",
        "точка": "Точ",
        "планина": "Пл",
        "село": "с.",
        "височини": "вис",
        "младши": "Мл",
        "станция": "С-я",
        "проход": "Прох",
        "баща": "Бщ"
    },
    "classifications": {
        "шофиране": "Шоф",
        "плавен": "Пл",
        "място": "Мя",
        "тераса": "Тер",
        "магистрала": "М-ла",
        "площад": "Пл",
        "пеш": "Пеш",
        "залив": "З-в",
        "пътека": "П-ка",
        "платно": "Пл",
        "улица": "Ул",
        "алея": "Ал",
        "пешеходна": "Пеш",
        "точка": "Тч",
        "задминаване": "Задм",
        "кръгово": "Кр",
        "връх": "Вр",
        "съд": "Сд",
        "булевард": "Бул",
        "път": "Път",
        "скоростна": "Скор",
        "мост": "Мо"
    },
    "directions": {
        "северозапад": "СЗ",
        "североизток": "СИ",
        "югозапад": "ЮЗ",
        "югоизток": "ЮИ",
        "север": "С",
        "изток": "И",
        "юг": "Ю"
    }
}

},{}],6:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "comunicacions": "Com.",
        "entitat de població": "Nucli",
        "disseminat": "Diss.",
        "cap de municipi": "Cap",
        "indret": "Indr.",
        "comarca": "Cca.",
        "relleu del litoral": "Lit.",
        "municipi": "Mun.",
        "xarxa hidrogràfica": "Curs Fluv.",
        "equipament": "Equip.",
        "orografia": "Orogr.",
        "barri": "Barri",
        "edificació": "Edif.",
        "edificació històrica": "Edif. Hist.",
        "entitat descentralitzada": "E.M.D.",
        "element hidrogràfic": "Hidr."
    },
    "classifications": {
        "rotonda": "Rot.",
        "carrerada": "Ca.",
        "jardí": "J.",
        "paratge": "Pge.",
        "pont": "Pont",
        "lloc": "Lloc",
        "rambla": "Rbla.",
        "cases": "Cses.",
        "barranc": "Bnc.",
        "plana": "Plana",
        "polígon": "Pol.",
        "muralla": "Mur.",
        "enllaç": "Ellaç",
        "antiga carretera": "Actra",
        "glorieta": "Glor.",
        "autovia": "Autv.",
        "prolongació": "Prol.",
        "calçada": "Cda.",
        "carretera": "Ctra.",
        "pujada": "Pda.",
        "torrent": "T.",
        "disseminat": "Disse",
        "barri": "B.",
        "cinturó": "Cinto",
        "passera": "Psera",
        "sender": "Send.",
        "carrer": "C.",
        "sèquia": "Sèq.",
        "blocs": "Bloc",
        "rambleta": "Rblt.",
        "partida": "Par.",
        "costa": "Cos.",
        "sector": "Sec.",
        "corraló": "Crral",
        "urbanització": "Urb.",
        "autopista": "Autp.",
        "grup": "Gr.",
        "platja": "Pja.",
        "jardins": "J.",
        "complex": "Comp.",
        "portals": "Ptals",
        "finca": "Fin.",
        "travessera": "Trav.",
        "plaça": "Pl.",
        "travessia": "Trv.",
        "polígon industrial": "PI.",
        "passatge": "Ptge.",
        "apartaments": "Apmt.",
        "mirador": "Mira.",
        "antic": "Antic",
        "accés": "Acc.",
        "colònia": "Col.",
        "corriol": "Crol.",
        "portal": "Ptal.",
        "porta": "Pta.",
        "port": "Port",
        "carreró": "Cró.",
        "riera": "Ra.",
        "circumval·lació": "Cval.",
        "baixada": "Bda.",
        "placeta": "Plta.",
        "escala": "Esc.",
        "gran via": "GV",
        "rial": "Rial",
        "conjunt": "Conj.",
        "avinguda": "Av.",
        "esplanada": "Esp.",
        "cantonada": "Cant.",
        "ronda": "Rda.",
        "corredor": "Cdor.",
        "drecera": "Drec.",
        "passadís": "Pdís.",
        "viaducte": "Vdct.",
        "passeig": "Pg.",
        "veïnat": "Veï."
    },
    "directions": {
        "sudest": "SE",
        "sudoest": "SO",
        "nordest": "NE",
        "nordoest": "NO",
        "est": "E",
        "nord": "N",
        "oest": "O",
        "sud": "S"
    }
}

},{}],7:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "skole": "Sk.",
        "ved": "v.",
        "centrum": "C.",
        "sankt": "Skt.",
        "vestre": "v.",
        "hospital": "Hosp.",
        "stræde": "Str.",
        "nordre": "Nr.",
        "plads": "Pl.",
        "universitet": "Uni.",
        "vænge": "vg.",
        "station": "St."
    },
    "classifications": {
        "avenue": "Ave",
        "gammel": "Gl.",
        "dronning": "Dronn.",
        "sønder": "Sdr.",
        "nørre": "Nr.",
        "vester": "V.",
        "vestre": "V.",
        "øster": "Ø.",
        "østre": "Ø.",
        "boulevard": "Boul."
    },
    "directions": {
        "sydøst": "SØ",
        "nordvest": "NV",
        "syd": "S",
        "nordøst": "NØ",
        "sydvest": "SV",
        "vest": "V",
        "nord": "N",
        "øst": "Ø"
    }
}

},{}],8:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {},
    "classifications": {},
    "directions": {
        "osten": "O",
        "nordosten": "NO",
        "süden": "S",
        "nordwest": "NW",
        "norden": "N",
        "südost": "SO",
        "südwest": "SW",
        "westen": "W"
    }
}

},{}],9:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "square": "Sq",
        "centre": "Ctr",
        "sister": "Sr",
        "lake": "Lk",
        "fort": "Ft",
        "route": "Rte",
        "william": "Wm",
        "national": "Nat’l",
        "junction": "Jct",
        "center": "Ctr",
        "saint": "St",
        "saints": "SS",
        "station": "Sta",
        "mount": "Mt",
        "junior": "Jr",
        "mountain": "Mtn",
        "heights": "Hts",
        "university": "Univ",
        "school": "Sch",
        "international": "Int’l",
        "apartments": "Apts",
        "crossing": "Xing",
        "creek": "Crk",
        "township": "Twp",
        "downtown": "Dtwn",
        "father": "Fr",
        "senior": "Sr",
        "point": "Pt",
        "river": "Riv",
        "market": "Mkt",
        "village": "Vil",
        "park": "Pk",
        "memorial": "Mem"
    },
    "classifications": {
        "place": "Pl",
        "circle": "Cir",
        "bypass": "Byp",
        "motorway": "Mwy",
        "crescent": "Cres",
        "road": "Rd",
        "cove": "Cv",
        "lane": "Ln",
        "square": "Sq",
        "street": "St",
        "freeway": "Fwy",
        "walk": "Wk",
        "plaza": "Plz",
        "parkway": "Pky",
        "avenue": "Ave",
        "pike": "Pk",
        "drive": "Dr",
        "highway": "Hwy",
        "footway": "Ftwy",
        "point": "Pt",
        "court": "Ct",
        "terrace": "Ter",
        "walkway": "Wky",
        "alley": "Aly",
        "expressway": "Expy",
        "bridge": "Br",
        "boulevard": "Blvd",
        "turnpike": "Tpk"
    },
    "directions": {
        "southeast": "SE",
        "northwest": "NW",
        "south": "S",
        "west": "W",
        "southwest": "SW",
        "north": "N",
        "east": "E",
        "northeast": "NE"
    }
}

},{}],10:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "segunda": "2ª",
        "octubre": "8bre",
        "doctores": "Drs",
        "doctora": "Dra",
        "internacional": "Intl",
        "doctor": "Dr",
        "segundo": "2º",
        "señorita": "Srta",
        "doctoras": "Drs",
        "primera": "1ª",
        "primero": "1º",
        "san": "S",
        "colonia": "Col",
        "doña": "Dña",
        "septiembre": "7bre",
        "diciembre": "10bre",
        "señor": "Sr",
        "ayuntamiento": "Ayto",
        "señora": "Sra",
        "tercera": "3ª",
        "tercero": "3º",
        "don": "D",
        "santa": "Sta",
        "ciudad": "Cdad",
        "noviembre": "9bre",
        "departamento": "Dep"
    },
    "classifications": {
        "camino": "Cmno",
        "avenida": "Av",
        "paseo": "Pº",
        "autopista": "Auto",
        "calle": "C",
        "plaza": "Pza",
        "carretera": "Crta"
    },
    "directions": {
        "este": "E",
        "noreste": "NE",
        "sur": "S",
        "suroeste": "SO",
        "noroeste": "NO",
        "oeste": "O",
        "sureste": "SE",
        "norte": "N"
    }
}

},{}],11:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "allée": "All",
        "aérodrome": "Aérod",
        "aéroport": "Aérop"
    },
    "classifications": {
        "centrale": "Ctrale",
        "campings": "Camp.",
        "urbains": "Urb.",
        "mineure": "Min.",
        "publique": "Publ.",
        "supérieur": "Sup.",
        "fédération": "Féd.",
        "notre-dame": "ND",
        "saint": "St",
        "centre hospitalier régional": "CHR",
        "exploitation": "Exploit.",
        "général": "Gal",
        "civiles": "Civ.",
        "maritimes": "Marit.",
        "aviation": "Aviat.",
        "iii": "3",
        "archéologique": "Archéo.",
        "musical": "Music.",
        "musicale": "Music.",
        "immeuble": "Imm.",
        "xv": "15",
        "hôtel": "Hôt.",
        "alpine": "Alp.",
        "communale": "Commun.",
        "v": "5",
        "global": "Glob.",
        "université": "Univ.",
        "confédéral": "Conféd.",
        "xx": "20",
        "x": "10",
        "piscine": "Pisc.",
        "dimanche": "di.",
        "fleuve": "Flv",
        "postaux": "Post.",
        "musicienne": "Music.",
        "département": "Dépt",
        "février": "Févr.",
        "municipales": "Munic.",
        "province": "Prov.",
        "communautés": "Commtés",
        "barrage": "Barr.",
        "mercredi": "me.",
        "présidentes": "Pdtes",
        "cafétérias": "Cafét.",
        "théâtral": "Thé.",
        "viticulteur": "Vitic.",
        "poste": "Post.",
        "spécialisée": "Spéc.",
        "agriculture": "Agric.",
        "infirmier": "Infirm.",
        "animation": "Anim.",
        "mondiale": "Mond.",
        "arrêt": "Arr.",
        "zone": "zon.",
        "municipaux": "Munic.",
        "grand": "Gd",
        "janvier": "Janv.",
        "fondateur": "Fond.",
        "première": "1re",
        "municipale": "Munic.",
        "direction": "Dir.",
        "anonyme": "Anon.",
        "départementale": "Dépt",
        "moyens": "Moy.",
        "novembre": "Nov.",
        "jardin": "Jard.",
        "petites": "Pet.",
        "privé": "Priv.",
        "centres": "Ctres",
        "forestier": "Forest.",
        "xiv": "14",
        "africaines": "Afric.",
        "sergent": "Sgt",
        "européenne": "Eur.",
        "privée": "Priv.",
        "café": "Cfé",
        "xix": "19",
        "hautes": "Htes",
        "major": "Mjr",
        "vendredi": "ve.",
        "municipalité": "Munic.",
        "sous-préfecture": "Ss-préf.",
        "spéciales": "Spéc.",
        "secondaires": "Second.",
        "viie": "7e",
        "moyenne": "Moy.",
        "commerciale": "Commerc.",
        "région": "Rég.",
        "américaines": "Amér.",
        "américains": "Amér.",
        "service": "Sce",
        "professeur": "Prof.",
        "départemental": "Dépt",
        "hôtels": "Hôt.",
        "mondiales": "Mond.",
        "ire": "1re",
        "caporal": "Capo.",
        "militaire": "Milit.",
        "lycée d'enseignement professionnel": "LEP",
        "adjudant": "Adj.",
        "médicale": "Méd.",
        "conférences": "Confér.",
        "universelle": "Univ.",
        "xiie": "12e",
        "supérieures": "Sup.",
        "naturel": "Natur.",
        "société nationale": "SN",
        "hospitalier": "Hosp.",
        "culturelle": "Cult.",
        "américain": "Amér.",
        "son altesse royale": "S.A.R.",
        "infirmière": "Infirm.",
        "viii": "8",
        "fondatrice": "Fond.",
        "madame": "Mme",
        "métropolitain": "Métrop.",
        "ophtalmologues": "Ophtalmos",
        "xviie": "18e",
        "viiie": "8e",
        "commerçante": "Commerç.",
        "centre d'enseignement du second degré": "CES",
        "septembre": "Sept.",
        "agriculteur": "Agric.",
        "xiii": "13",
        "pontifical": "Pontif.",
        "cafétéria": "Cafét.",
        "prince": "Pce",
        "vie": "6e",
        "archiduchesse": "Archid.",
        "occidental": "Occ.",
        "spectacles": "Spect.",
        "camping": "Camp.",
        "métro": "Mº",
        "arrondissement": "Arrond.",
        "viticole": "Vitic.",
        "ii": "2",
        "siècle": "Si.",
        "chapelles": "Chap.",
        "centre": "Ctre",
        "sapeur-pompiers": "Sap.-pomp.",
        "établissements": "Étabts",
        "société anonyme": "SA",
        "directeurs": "Dir.",
        "vii": "7",
        "culturel": "Cult.",
        "central": "Ctral",
        "métropolitaine": "Métrop.",
        "administrations": "Admin.",
        "amiraux": "Amir.",
        "sur": "s/",
        "premiers": "1ers",
        "provence-alpes-côte d'azur": "PACA",
        "cathédrale": "Cathéd.",
        "iv": "4",
        "postale": "Post.",
        "social": "Soc.",
        "spécialisé": "Spéc.",
        "district": "Distr.",
        "technologique": "Techno.",
        "viticoles": "Vitic.",
        "ix": "9",
        "protégés": "Prot.",
        "historiques": "Hist.",
        "sous": "s/s",
        "national": "Nal",
        "ambassade": "Amb.",
        "cafés": "Cfés",
        "agronomie": "Agro.",
        "sapeurs": "Sap.",
        "petits": "Pet.",
        "monsieur": "M.",
        "boucher": "Bouch.",
        "restaurant": "Restau.",
        "lycée": "Lyc.",
        "urbaine": "Urb.",
        "préfecture": "Préf.",
        "districts": "Distr.",
        "civil": "Civ.",
        "protégées": "Prot.",
        "sapeur": "Sap.",
        "théâtre": "Thé.",
        "collège": "Coll.",
        "mardi": "ma.",
        "mémorial": "Mémor.",
        "africain": "Afric.",
        "républicaine": "Républ.",
        "sociale": "Soc.",
        "spécial": "Spéc.",
        "technologie": "Techno.",
        "charcuterie": "Charc.",
        "commerces": "Commerc.",
        "fluviale": "Flv",
        "parachutistes": "Para.",
        "primaires": "Prim.",
        "directions": "Dir.",
        "présidentiel": "Pdtl",
        "nationales": "Nales",
        "après": "apr.",
        "samedi": "sa.",
        "unité": "U.",
        "xxiii": "23",
        "associé": "Assoc.",
        "électrique": "Électr.",
        "populaire": "Pop.",
        "asiatique": "Asiat.",
        "navigable": "Navig.",
        "présidente": "Pdte",
        "xive": "14e",
        "associés": "Assoc.",
        "pompiers": "Pomp.",
        "agricoles": "Agric.",
        "élém": "Élém.",
        "décembre": "Déc.",
        "son altesse": "S.Alt.",
        "après-midi": "a.-m.",
        "mineures": "Min.",
        "juillet": "Juil.",
        "aviatrices": "Aviat.",
        "fondation": "Fond.",
        "pontificaux": "Pontif.",
        "temple": "Tple",
        "européennes": "Eur.",
        "régionale": "Rég.",
        "informations": "Infos",
        "mondiaux": "Mond.",
        "infanterie": "Infant.",
        "archéologie": "Archéo.",
        "dans": "d/",
        "hospice": "Hosp.",
        "spectacle": "Spect.",
        "hôtels-restaurants": "Hôt.-Rest.",
        "hôtel-restaurant": "Hôt.-Rest.",
        "hélicoptère": "hélico",
        "xixe": "19e",
        "cliniques": "Clin.",
        "docteur": "Dr",
        "secondaire": "Second.",
        "municipal": "Munic.",
        "générale": "Gale",
        "château": "Chât.",
        "commerçant": "Commerç.",
        "avril": "Avr.",
        "clinique": "Clin.",
        "urbaines": "Urb.",
        "navale": "Nav.",
        "navigation": "Navig.",
        "asiatiques": "Asiat.",
        "pontificales": "Pontif.",
        "administrative": "Admin.",
        "syndicat": "Synd.",
        "lundi": "lu.",
        "petite": "Pet.",
        "maritime": "Marit.",
        "métros": "Mº",
        "enseignement": "Enseign.",
        "fluviales": "Flv",
        "historique": "Hist.",
        "comtés": "Ctés",
        "résidentiel": "Résid.",
        "international": "Int.",
        "supérieure": "Sup.",
        "centre hospitalier universitaire": "CHU",
        "confédération": "Conféd.",
        "boucherie": "Bouch.",
        "fondatrices": "Fond.",
        "médicaux": "Méd.",
        "européens": "Eur.",
        "orientaux": "Ori.",
        "naval": "Nav.",
        "étang": "Étg",
        "provincial": "Prov.",
        "junior": "Jr",
        "départementales": "Dépt",
        "musique": "Musiq.",
        "directrices": "Dir.",
        "maréchal": "Mal",
        "civils": "Civ.",
        "protégé": "Prot.",
        "établissement": "Étabt",
        "trafic": "Traf.",
        "aviateur": "Aviat.",
        "archives": "Arch.",
        "africains": "Afric.",
        "maternelle": "Matern.",
        "industrielle": "Ind.",
        "administratif": "Admin.",
        "oriental": "Ori.",
        "universitaire": "Univ.",
        "majeur": "Maj.",
        "haute": "Hte",
        "communal": "Commun.",
        "petit": "Pet.",
        "commune": "Commun.",
        "exploitant": "Exploit.",
        "conférence": "Confér.",
        "monseigneur": "Mgr",
        "pharmacien": "Pharm.",
        "jeudi": "je.",
        "primaire": "Prim.",
        "hélicoptères": "hélicos",
        "agronomique": "Agro.",
        "médecin": "Méd.",
        "ve": "5e",
        "pontificale": "Pontif.",
        "ier": "1er",
        "cinéma": "Ciné",
        "fluvial": "Flv",
        "occidentaux": "Occ.",
        "commerçants": "Commerç.",
        "banque": "Bq",
        "moyennes": "Moy.",
        "pharmacienne": "Pharm.",
        "démocratique": "Dém.",
        "cinémas": "Cinés",
        "spéciale": "Spéc.",
        "présidents": "Pdts",
        "directrice": "Dir.",
        "vi": "6",
        "basse": "Bas.",
        "xve": "15e",
        "état": "É.",
        "aviateurs": "Aviat.",
        "majeurs": "Maj.",
        "infirmiers": "Infirm.",
        "église": "Égl.",
        "confédérale": "Conféd.",
        "xxie": "21e",
        "comte": "Cte",
        "européen": "Eur.",
        "union": "U.",
        "pharmacie": "Pharm.",
        "infirmières": "Infirm.",
        "comté": "Cté",
        "sportive": "Sport.",
        "deuxième": "2e",
        "xvi": "17",
        "haut": "Ht",
        "médicales": "Méd.",
        "développé": "Dévelop.",
        "bâtiment": "Bât.",
        "commerce": "Commerc.",
        "ive": "4e",
        "associatif": "Assoc.",
        "rural": "Rur.",
        "cimetière": "Cim.",
        "régional": "Rég.",
        "ferroviaire": "Ferr.",
        "vers": "v/",
        "mosquée": "Mosq.",
        "mineurs": "Min.",
        "nautique": "Naut.",
        "châteaux": "Chât.",
        "sportif": "Sport.",
        "mademoiselle": "Mle",
        "école": "Éc.",
        "doyen": "Doy.",
        "industriel": "Ind.",
        "chapelle": "Chap.",
        "sociétés": "Stés",
        "internationale": "Int.",
        "coopératif": "Coop.",
        "hospices": "Hosp.",
        "xxii": "22",
        "parachutiste": "Para.",
        "alpines": "Alp.",
        "civile": "Civ.",
        "xvie": "17e",
        "états": "É.",
        "musée": "Msée",
        "centrales": "Ctrales",
        "globaux": "Glob.",
        "supérieurs": "Sup.",
        "syndicats": "Synd.",
        "archevêque": "Archev.",
        "docteurs": "Drs",
        "bibliothèque": "Biblio.",
        "lieutenant": "Lieut.",
        "république": "Rép.",
        "vétérinaire": "Vét.",
        "départementaux": "Dépt",
        "premier": "1er",
        "fluviaux": "Flv",
        "animé": "Anim.",
        "orientales": "Ori.",
        "technologiques": "Techno.",
        "princesse": "Pse",
        "routière": "Rout.",
        "coopérative": "Coop.",
        "scolaire": "Scol.",
        "écoles": "Éc.",
        "football": "Foot",
        "territoriale": "Territ.",
        "commercial": "Commerc.",
        "mineur": "Min.",
        "millénaires": "Mill.",
        "association": "Assoc.",
        "catholique": "Cathol.",
        "administration": "Admin.",
        "mairie": "Mair.",
        "portuaire": "Port.",
        "tertiaires": "Terti.",
        "théâtrale": "Thé.",
        "palais": "Pal.",
        "troisième": "3e",
        "directeur": "Dir.",
        "vétérinaires": "Vét.",
        "faculté": "Fac.",
        "occidentales": "Occ.",
        "viticulteurs": "Vitic.",
        "xvii": "18",
        "occidentale": "Occ.",
        "amiral": "Amir.",
        "professionnel": "Profess.",
        "administratives": "Admin.",
        "commerciales": "Commerc.",
        "saints": "Sts",
        "agronomes": "Agro.",
        "stade": "Std",
        "sous-préfet": "Ss-préf.",
        "senior": "Sr",
        "agronome": "Agro.",
        "terrain": "Terr.",
        "catholiques": "Cathol.",
        "résidentielle": "Résid.",
        "grands": "Gds",
        "exploitants": "Exploit.",
        "xiiie": "13e",
        "croix": "Cx",
        "généraux": "Gaux",
        "crédit": "Créd.",
        "cimetières": "Cim.",
        "antenne": "Ant.",
        "médical": "Méd.",
        "collèges": "Coll.",
        "musicien": "Music.",
        "apostolique": "Apost.",
        "postal": "Post.",
        "territorial": "Territ.",
        "urbanisme": "Urb.",
        "préfectorale": "Préf.",
        "fondateurs": "Fond.",
        "information": "Info.",
        "églises": "Égl.",
        "ophtalmologue": "Ophtalmo",
        "congrégation": "Congrég.",
        "charcutier": "Charc.",
        "étage": "ét.",
        "consulat": "Consul.",
        "public": "Publ.",
        "ferrée": "Ferr.",
        "matin": "mat.",
        "société anonyme à responsabilité limitée": "SARL",
        "monuments": "Mmts",
        "protection": "Prot.",
        "universel": "Univ.",
        "nationale": "Nale",
        "président": "Pdt",
        "provinciale": "Prov.",
        "agriculteurs": "Agric.",
        "préfectoral": "Préf.",
        "xxe": "20e",
        "alpins": "Alp.",
        "avant": "av.",
        "infirmerie": "Infirm.",
        "deux mil": "2000",
        "rurale": "Rur.",
        "administratifs": "Admin.",
        "octobre": "Oct.",
        "archipel": "Archip.",
        "communauté": "Commté",
        "globales": "Glob.",
        "alpin": "Alp.",
        "numéros": "Nºˢ",
        "lieutenant-colonel": "Lieut.-Col.",
        "jésus-christ": "J.-C.",
        "agricole": "Agric.",
        "sa majesté": "S.Maj.",
        "associative": "Assoc.",
        "xxi": "21",
        "présidentielle": "Pdtle",
        "moyen": "Moy.",
        "fédéral": "Féd.",
        "professionnelle": "Profess.",
        "tertiaire": "Terti.",
        "ixe": "9e",
        "hôpital": "Hôp.",
        "technologies": "Techno.",
        "iiie": "3e",
        "développement": "Dévelop.",
        "monument": "Mmt",
        "forestière": "Forest.",
        "numéro": "Nº",
        "viticulture": "Vitic.",
        "traversière": "Traver.",
        "technique": "Tech.",
        "électriques": "Électr.",
        "militaires": "Milit.",
        "pompier": "Pomp.",
        "américaine": "Amér.",
        "préfet": "Préf.",
        "congrégations": "Congrég.",
        "pâtissier": "Pâtiss.",
        "mondial": "Mond.",
        "ophtalmologie": "Ophtalm.",
        "sainte": "Ste",
        "africaine": "Afric.",
        "aviatrice": "Aviat.",
        "doyens": "Doy.",
        "société": "Sté",
        "majeures": "Maj.",
        "orientale": "Ori.",
        "ministère": "Min.",
        "archiduc": "Archid.",
        "territoire": "Territ.",
        "techniques": "Tech.",
        "île-de-france": "IDF",
        "globale": "Glob.",
        "xe": "10e",
        "xie": "11e",
        "majeure": "Maj.",
        "commerciaux": "Commerc.",
        "maire": "Mair.",
        "spéciaux": "Spéc.",
        "grande": "Gde",
        "messieurs": "MM",
        "colonel": "Col.",
        "millénaire": "Mill.",
        "xi": "11",
        "urbain": "Urb.",
        "fédérale": "Féd.",
        "ferré": "Ferr.",
        "rivière": "Riv.",
        "républicain": "Républ.",
        "grandes": "Gdes",
        "régiment": "Régim.",
        "hauts": "Hts",
        "catégorie": "Catég.",
        "basses": "Bas.",
        "xii": "12",
        "agronomiques": "Agro.",
        "iie": "2e",
        "protégée": "Prot.",
        "sapeur-pompier": "Sap.-pomp."
    },
    "directions": {
        "est-nord-est": "ENE",
        "nord-est": "NE",
        "ouest": "O",
        "sud-est": "SE",
        "est-sud-est": "ESE",
        "nord-nord-est": "NNE",
        "sud": "S",
        "nord-nord-ouest": "NNO",
        "nord-ouest": "NO",
        "nord": "N",
        "ouest-sud-ouest": "OSO",
        "ouest-nord-ouest": "ONO",
        "sud-ouest": "SO",
        "sud-sud-est": "SSE",
        "sud-sud-ouest": "SSO",
        "est": "E"
    }
}

},{}],12:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "שדרות": "שד'"
    },
    "classifications": {},
    "directions": {}
}

},{}],13:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {},
    "classifications": {},
    "directions": {
        "kelet": "K",
        "északkelet": "ÉK",
        "dél": "D",
        "északnyugat": "ÉNY",
        "észak": "É",
        "délkelet": "DK",
        "délnyugat": "DNY",
        "nyugat": "NY"
    }
}

},{}],14:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "apartamentai": "Apt",
        "aukštumos": "Aukš",
        "centras": "Ctr",
        "ežeras": "Ež",
        "fortas": "Ft",
        "greitkelis": "Grtkl",
        "juosta": "Jst",
        "kaimas": "Km",
        "kalnas": "Kln",
        "kelias": "Kl",
        "kiemelis": "Kml",
        "miestelis": "Mstl",
        "miesto centras": "M.Ctr",
        "mokykla": "Mok",
        "nacionalinis": "Nac",
        "paminklas": "Pmkl",
        "parkas": "Pk",
        "pusratis": "Psrt",
        "sankryža": "Skrž",
        "sesė": "Sesė",
        "skveras": "Skv",
        "stotis": "St",
        "šv": "Šv",
        "tarptautinis": "Trptaut",
        "taškas": "Tšk",
        "tėvas": "Tėv",
        "turgus": "Tgs",
        "universitetas": "Univ",
        "upė": "Up",
        "upelis": "Up",
        "vieta": "Vt"
    },
    "classifications": {
        "aikštė": "a.",
        "alėja": "al.",
        "aplinkkelis": "aplinkl.",
        "autostrada": "auto.",
        "bulvaras": "b.",
        "gatvė": "g.",
        "kelias": "kel.",
        "krantinė": "krant.",
        "prospektas": "pr.",
        "plentas": "pl.",
        "skersgatvis": "skg.",
        "takas": "tak.",
        "tiltas": "tlt."
    },
    "directions": {
        "pietūs": "P",
        "vakarai": "V",
        "šiaurė": "Š",
        "šiaurės vakarai": "ŠV",
        "pietryčiai": "PR",
        "šiaurės rytai": "ŠR",
        "rytai": "R",
        "pietvakariai": "PV"
    }
}

},{}],15:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "centrum": "Cntrm",
        "nationaal": "Nat’l",
        "berg": "Brg",
        "meer": "Mr",
        "kruising": "Krsng",
        "toetreden": "Ttrdn"
    },
    "classifications": {
        "bypass": "Pass",
        "brug": "Br",
        "straat": "Str",
        "rechtbank": "Rbank",
        "snoek": "Snk",
        "autobaan": "Baan",
        "terras": "Trrs",
        "punt": "Pt",
        "plaza": "Plz",
        "rijden": "Rijd",
        "parkway": "Pky",
        "inham": "Nham",
        "snelweg": "Weg",
        "halve maan": "Maan",
        "cirkel": "Crkl",
        "laan": "Ln",
        "rijbaan": "Strook",
        "weg": "Weg",
        "lopen": "Lpn",
        "autoweg": "Weg",
        "boulevard": "Blvd",
        "plaats": "Plts",
        "steeg": "Stg",
        "voetpad": "Stoep"
    },
    "directions": {
        "noordoost": "NO",
        "westen": "W",
        "zuiden": "Z",
        "zuidwest": "ZW",
        "oost": "O",
        "zuidoost": "ZO",
        "noordwest": "NW",
        "noorden": "N"
    }
}

},{}],16:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "апостола": "ап.",
        "апостолов": "апп.",
        "великомученика": "вмч",
        "великомученицы": "вмц.",
        "владение": "вл.",
        "город": "г.",
        "деревня": "д.",
        "имени": "им.",
        "мученика":"мч.",
        "мучеников": "мчч.",
        "мучениц": "мцц.",
        "мученицы": "мц.",
        "озеро": "о.",
        "посёлок": "п.",
        "преподобного":  "прп.",
        "преподобных": "прпп.",
        "река": "р.",
        "святителей": "свтт.",
        "святителя": "свт.",
        "священномученика": "сщмч.",
        "священномучеников": "сщмчч.",
        "станция": "ст.",
        "участок": "уч."
    },
    "classifications": {
        "проезд": "пр-д",
        "проспект": "пр.",
        "переулок": "пер.",
        "набережная": "наб.",
        "площадь": "пл.",
        "шоссе": "ш.",
        "бульвар": "б.",
        "тупик": "туп.",
        "улица": "ул."
    },
    "directions": {
        "восток": "В",
        "северо-восток": "СВ",
        "юго-восток": "ЮВ",
        "юго-запад": "ЮЗ",
        "северо-запад": "СЗ",
        "север": "С",
        "запад": "З",
        "юг": "Ю"
    }
}

},{}],17:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {},
    "classifications": {},
    "directions": {
        "vzhod": "V",
        "severovzhod": "SV",
        "jug": "J",
        "severozahod": "SZ",
        "sever": "S",
        "jugovzhod": "JV",
        "jugozahod": "JZ",
        "zahod": "Z"
    }
}

},{}],18:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "sankta": "s:ta",
        "gamla": "G:la",
        "sankt": "s:t"
    },
    "classifications": {
        "Bro": "Br"
    },
    "directions": {
        "norr": "N",
        "sydöst": "SO",
        "väster": "V",
        "öster": "O",
        "nordväst": "NV",
        "sydväst": "SV",
        "söder": "S",
        "nordöst": "NO"
    }
}

},{}],19:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {},
    "classifications": {},
    "directions": {
        "схід": "Сх",
        "північний схід": "ПнСх",
        "південь": "Пд",
        "північний захід": "ПнЗд",
        "північ": "Пн",
        "південний схід": "ПдСх",
        "південний захід": "ПдЗх",
        "захід": "Зх"
    }
}

},{}],20:[function(_dereq_,module,exports){
module.exports={
    "abbreviations": {
        "viện bảo tàng": "VBT",
        "thị trấn": "Tt",
        "đại học": "ĐH",
        "căn cứ không quan": "CCKQ",
        "câu lạc bộ": "CLB",
        "bưu điện": "BĐ",
        "khách sạn": "KS",
        "khu du lịch": "KDL",
        "khu công nghiệp": "KCN",
        "khu nghỉ mát": "KNM",
        "thị xã": "Tx",
        "khu chung cư": "KCC",
        "phi trường": "PT",
        "trung tâm": "TT",
        "tổng công ty": "TCty",
        "trung học cơ sở": "THCS",
        "sân bay quốc tế": "SBQT",
        "trung học phổ thông": "THPT",
        "cao đẳng": "CĐ",
        "công ty": "Cty",
        "sân bay": "SB",
        "thành phố": "Tp",
        "công viên": "CV",
        "sân vận động": "SVĐ",
        "linh mục": "LM",
        "vườn quốc gia": "VQG"
    },
    "classifications": {
        "huyện lộ": "HL",
        "đường tỉnh": "ĐT",
        "quốc lộ": "QL",
        "xa lộ": "XL",
        "hương lộ": "HL",
        "tỉnh lộ": "TL",
        "đường huyện": "ĐH",
        "đường cao tốc": "ĐCT",
        "đại lộ": "ĐL",
        "việt nam": "VN",
        "quảng trường": "QT",
        "đường bộ": "ĐB"
    },
    "directions": {
        "tây": "T",
        "nam": "N",
        "đông nam": "ĐN",
        "đông bắc": "ĐB",
        "tây nam": "TN",
        "đông": "Đ",
        "bắc": "B"
    }
}

},{}],21:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "regExpFlags": "gi"
    },
    "v5": {
        "article": [
            ["^ Acc[èe]s ", " l’accès "],
            ["^ Aire ", " l’aire "],
            ["^ All[ée]e ", " l’allée "],
            ["^ Anse ", " l’anse "],
            ["^ (L['’])?Autoroute ", " l’autoroute "],
            ["^ Avenue ", " l’avenue "],
            ["^ Barreau ", " le barreau "],
            ["^ Boulevard ", " le boulevard "],
            ["^ Chemin ", " le chemin "],
            ["^ Petit[\\- ]Chemin ", " le petit chemin "],
            ["^ Cit[ée] ", " la cité "],
            ["^ Clos ", " le clos "],
            ["^ Corniche ", " la corniche "],
            ["^ Cour ", " la cour "],
            ["^ Cours ", " le cours "],
            ["^ D[ée]viation ", " la déviation "],
            ["^ Entr[ée]e ", " l’entrée "],
            ["^ Esplanade ", " l’esplanade "],
            ["^ Galerie ", " la galerie "],
            ["^ Impasse ", " l’impasse "],
            ["^ Lotissement ", " le lotissement "],
            ["^ Mont[ée]e ", " la montée "],
            ["^ Parc ", " le parc "],
            ["^ Parvis ", " le parvis "],
            ["^ Passage ", " le passage "],
            ["^ Place ", " la place "],
            ["^ Petit[\\- ]Pont ", " le petit-pont "],
            ["^ Pont ", " le pont "],
            ["^ Promenade ", " la promenade "],
            ["^ Quai ", " le quai "],
            ["^ Rocade ", " la rocade "],
            ["^ Rond[\\- ]?Point ", " le rond-point "],
            ["^ Route ", " la route "],
            ["^ Rue ", " la rue "],
            ["^ Grande Rue ", " la grande rue "],
            ["^ Sente ", " la sente "],
            ["^ Sentier ", " le sentier "],
            ["^ Sortie ", " la sortie "],
            ["^ Souterrain ", " le souterrain "],
            ["^ Square ", " le square "],
            ["^ Terrasse ", " la terrasse "],
            ["^ Traverse ", " la traverse "],
            ["^ Tunnel ", " le tunnel "],
            ["^ Viaduc ", " le viaduc "],
            ["^ Villa ", " la villa "],
            ["^ Village ", " le village "],
            ["^ Voie ", " la voie "],

            [" ([dl])'", " $1’"]
        ],
        "preposition": [
            ["^ Le ", "  du "],
            ["^ Les ", "  des "],
            ["^ La ", "  de La "],

            ["^ Acc[èe]s ", "  de l’accès "],
            ["^ Aire ", "  de l’aire "],
            ["^ All[ée]e ", "  de l’allée "],
            ["^ Anse ", "  de l’anse "],
            ["^ (L['’])?Autoroute ", "  de l’autoroute "],
            ["^ Avenue ", "  de l’avenue "],
            ["^ Barreau ", "  du barreau "],
            ["^ Boulevard ", "  du boulevard "],
            ["^ Chemin ", "  du chemin "],
            ["^ Petit[\\- ]Chemin ", "  du petit chemin "],
            ["^ Cit[ée] ", "  de la cité "],
            ["^ Clos ", "  du clos "],
            ["^ Corniche ", "  de la corniche "],
            ["^ Cour ", "  de la cour "],
            ["^ Cours ", "  du cours "],
            ["^ D[ée]viation ", "  de la déviation "],
            ["^ Entr[ée]e ", "  de l’entrée "],
            ["^ Esplanade ", "  de l’esplanade "],
            ["^ Galerie ", "  de la galerie "],
            ["^ Impasse ", "  de l’impasse "],
            ["^ Lotissement ", "  du lotissement "],
            ["^ Mont[ée]e ", "  de la montée "],
            ["^ Parc ", "  du parc "],
            ["^ Parvis ", "  du parvis "],
            ["^ Passage ", "  du passage "],
            ["^ Place ", "  de la place "],
            ["^ Petit[\\- ]Pont ", "  du petit-pont "],
            ["^ Pont ", "  du pont "],
            ["^ Promenade ", "  de la promenade "],
            ["^ Quai ", "  du quai "],
            ["^ Rocade ", "  de la rocade "],
            ["^ Rond[\\- ]?Point ", "  du rond-point "],
            ["^ Route ", "  de la route "],
            ["^ Rue ", "  de la rue "],
            ["^ Grande Rue ", "  de la grande rue "],
            ["^ Sente ", "  de la sente "],
            ["^ Sentier ", "  du sentier "],
            ["^ Sortie ", "  de la sortie "],
            ["^ Souterrain ", "  du souterrain "],
            ["^ Square ", "  du square "],
            ["^ Terrasse ", "  de la terrasse "],
            ["^ Traverse ", "  de la traverse "],
            ["^ Tunnel ", "  du tunnel "],
            ["^ Viaduc ", "  du viaduc "],
            ["^ Villa ", "  de la villa "],
            ["^ Village ", "  du village "],
            ["^ Voie ", "  de la voie "],

            ["^ ([AÂÀEÈÉÊËIÎÏOÔUÙÛÜYŸÆŒ])", "  d’$1"],
            ["^ (\\S)", "  de $1"],
            [" ([dl])'", " $1’"]
        ],
        "rotary": [
            ["^ Le ", "  le rond-point du "],
            ["^ Les ", "  le rond-point des "],
            ["^ La ", "  le rond-point de La "],

            ["^ Acc[èe]s ", " le rond-point de l’accès "],
            ["^ Aire ", "  le rond-point de l’aire "],
            ["^ All[ée]e ", "  le rond-point de l’allée "],
            ["^ Anse ", "  le rond-point de l’anse "],
            ["^ (L['’])?Autoroute ", "  le rond-point de l’autoroute "],
            ["^ Avenue ", "  le rond-point de l’avenue "],
            ["^ Barreau ", "  le rond-point du barreau "],
            ["^ Boulevard ", "  le rond-point du boulevard "],
            ["^ Chemin ", "  le rond-point du chemin "],
            ["^ Petit[\\- ]Chemin ", "  le rond-point du petit chemin "],
            ["^ Cit[ée] ", "  le rond-point de la cité "],
            ["^ Clos ", "  le rond-point du clos "],
            ["^ Corniche ", "  le rond-point de la corniche "],
            ["^ Cour ", "  le rond-point de la cour "],
            ["^ Cours ", "  le rond-point du cours "],
            ["^ D[ée]viation ", "  le rond-point de la déviation "],
            ["^ Entr[ée]e ", "  le rond-point de l’entrée "],
            ["^ Esplanade ", "  le rond-point de l’esplanade "],
            ["^ Galerie ", "  le rond-point de la galerie "],
            ["^ Impasse ", "  le rond-point de l’impasse "],
            ["^ Lotissement ", "  le rond-point du lotissement "],
            ["^ Mont[ée]e ", "  le rond-point de la montée "],
            ["^ Parc ", "  le rond-point du parc "],
            ["^ Parvis ", "  le rond-point du parvis "],
            ["^ Passage ", "  le rond-point du passage "],
            ["^ Place ", "  le rond-point de la place "],
            ["^ Petit[\\- ]Pont ", "  le rond-point du petit-pont "],
            ["^ Pont ", "  le rond-point du pont "],
            ["^ Promenade ", "  le rond-point de la promenade "],
            ["^ Quai ", "  le rond-point du quai "],
            ["^ Rocade ", "  le rond-point de la rocade "],
            ["^ Rond[\\- ]?Point ", "  le rond-point "],
            ["^ Route ", "  le rond-point de la route "],
            ["^ Rue ", "  le rond-point de la rue "],
            ["^ Grande Rue ", "  le rond-point de la grande rue "],
            ["^ Sente ", "  le rond-point de la sente "],
            ["^ Sentier ", "  le rond-point du sentier "],
            ["^ Sortie ", "  le rond-point de la sortie "],
            ["^ Souterrain ", "  le rond-point du souterrain "],
            ["^ Square ", "  le rond-point du square "],
            ["^ Terrasse ", "  le rond-point de la terrasse "],
            ["^ Traverse ", "  le rond-point de la traverse "],
            ["^ Tunnel ", "  le rond-point du tunnel "],
            ["^ Viaduc ", "  le rond-point du viaduc "],
            ["^ Villa ", "  le rond-point de la villa "],
            ["^ Village ", "  le rond-point du village "],
            ["^ Voie ", "  le rond-point de la voie "],

            ["^ ([AÂÀEÈÉÊËIÎÏOÔUÙÛÜYŸÆŒ])", "  le rond-point d’$1"],
            ["^ (\\S)", "  le rond-point de $1"],
            [" ([dl])'", " $1’"]
        ],
        "arrival": [
            ["^ Le ", "  au "],
            ["^ Les ", "  aux "],
            ["^ La ", "  à La "],
            ["^ (\\S)", "  à $1"],

            [" ([dl])'", " $1’"]
        ]
    }
}

},{}],22:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "regExpFlags": ""
    },
    "v5": {
        "accusative": [
            ["^ ([«\"])", " трасса $1"],

            ["^ (\\S+)ая [Аа]ллея ", " $1ую аллею "],
            ["^ (\\S+)ья [Аа]ллея ", " $1ью аллею "],
            ["^ (\\S+)яя [Аа]ллея ", " $1юю аллею "],
            ["^ (\\d+)-я (\\S+)ая [Аа]ллея ", " $1-ю $2ую аллею "],
            ["^ [Аа]ллея ", " аллею "],

            ["^ (\\S+)ая-(\\S+)ая [Уу]лица ", " $1ую-$2ую улицу "],
            ["^ (\\S+)ая [Уу]лица ", " $1ую улицу "],
            ["^ (\\S+)ья [Уу]лица ", " $1ью улицу "],
            ["^ (\\S+)яя [Уу]лица ", " $1юю улицу "],
            ["^ (\\d+)-я [Уу]лица ", " $1-ю улицу "],
            ["^ (\\d+)-я (\\S+)ая [Уу]лица ", " $1-ю $2ую улицу "],
            ["^ (\\S+)ая (\\S+)ая [Уу]лица ", " $1ую $2ую улицу "],
            ["^ (\\S+[вн])а [Уу]лица ", " $1у улицу "],
            ["^ (\\S+)ая (\\S+[вн])а [Уу]лица ", " $1ую $2у улицу "],
            ["^ Даньславля [Уу]лица ", " Даньславлю улицу "],
            ["^ Добрыня [Уу]лица ", " Добрыню улицу "],
            ["^ Людогоща [Уу]лица ", " Людогощу улицу "],
            ["^ [Уу]лица ", " улицу "],

            ["^ (\\d+)-я [Лл]иния ", " $1-ю линию "],
            ["^ (\\d+)-(\\d+)-я [Лл]иния ", " $1-$2-ю линию "],
            ["^ (\\S+)ая [Лл]иния ", " $1ую линию "],
            ["^ (\\S+)ья [Лл]иния ", " $1ью линию "],
            ["^ (\\S+)яя [Лл]иния ", " $1юю линию "],
            ["^ (\\d+)-я (\\S+)ая [Лл]иния ", " $1-ю $2ую линию "],
            ["^ [Лл]иния ", " линию "],

            ["^ (\\d+)-(\\d+)-я [Лл]инии ", " $1-$2-ю линии "],

            ["^ (\\S+)ая [Нн]абережная ", " $1ую набережную "],
            ["^ (\\S+)ья [Нн]абережная ", " $1ью набережную "],
            ["^ (\\S+)яя [Нн]абережная ", " $1юю набережную "],
            ["^ (\\d+)-я (\\S+)ая [Нн]абережная ", " $1-ю $2ую набережную "],
            ["^ [Нн]абережная ", " набережную "],

            ["^ (\\S+)ая [Пп]лощадь ", " $1ую площадь "],
            ["^ (\\S+)ья [Пп]лощадь ", " $1ью площадь "],
            ["^ (\\S+)яя [Пп]лощадь ", " $1юю площадь "],
            ["^ (\\S+[вн])а [Пп]лощадь ", " $1у площадь "],
            ["^ (\\d+)-я (\\S+)ая [Пп]лощадь ", " $1-ю $2ую площадь "],
            ["^ [Пп]лощадь ", " площадь "],

            ["^ (\\S+)ая [Пп]росека ", " $1ую просеку "],
            ["^ (\\S+)ья [Пп]росека ", " $1ью просеку "],
            ["^ (\\S+)яя [Пп]росека ", " $1юю просеку "],
            ["^ (\\d+)-я [Пп]росека ", " $1-ю просеку "],
            ["^ [Пп]росека ", " просеку "],

            ["^ (\\S+)ая [Ээ]стакада ", " $1ую эстакаду "],
            ["^ (\\S+)ья [Ээ]стакада ", " $1ью эстакаду "],
            ["^ (\\S+)яя [Ээ]стакада ", " $1юю эстакаду "],
            ["^ (\\d+)-я (\\S+)ая [Ээ]стакада ", " $1-ю $2ую эстакаду "],
            ["^ [Ээ]стакада ", " эстакаду "],

            ["^ (\\S+)ая [Мм]агистраль ", " $1ую магистраль "],
            ["^ (\\S+)ья [Мм]агистраль ", " $1ью магистраль "],
            ["^ (\\S+)яя [Мм]агистраль ", " $1юю магистраль "],
            ["^ (\\S+)ая (\\S+)ая [Мм]агистраль ", " $1ую $2ую магистраль "],
            ["^ (\\d+)-я (\\S+)ая [Мм]агистраль ", " $1-ю $2ую магистраль "],
            ["^ [Мм]агистраль ", " магистраль "],

            ["^ (\\S+)ая [Рр]азвязка ", " $1ую развязку "],
            ["^ (\\S+)ья [Рр]азвязка ", " $1ью развязку "],
            ["^ (\\S+)яя [Рр]азвязка ", " $1юю развязку "],
            ["^ (\\d+)-я (\\S+)ая [Рр]азвязка ", " $1-ю $2ую развязку "],
            ["^ [Рр]азвязка ", " развязку "],

            ["^ (\\S+)ая [Тт]расса ", " $1ую трассу "],
            ["^ (\\S+)ья [Тт]расса ", " $1ью трассу "],
            ["^ (\\S+)яя [Тт]расса ", " $1юю трассу "],
            ["^ (\\d+)-я (\\S+)ая [Тт]расса ", " $1-ю $2ую трассу "],
            ["^ [Тт]расса ", " трассу "],

            ["^ (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ую $2дорогу "],
            ["^ (\\S+)ья ([Аа]вто)?[Дд]орога ", " $1ью $2дорогу "],
            ["^ (\\S+)яя ([Аа]вто)?[Дд]орога ", " $1юю $2дорогу "],
            ["^ (\\S+)ая (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ую $2ую $3дорогу "],
            ["^ (\\d+)-я (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1-ю $2ую $3дорогу "],
            ["^ ([Аа]вто)?[Дд]орога ", " $1дорогу "],

            ["^ (\\S+)ая [Дд]орожка ", " $1ую дорожку "],
            ["^ (\\S+)ья [Дд]орожка ", " $1ью дорожку "],
            ["^ (\\S+)яя [Дд]орожка ", " $1юю дорожку "],
            ["^ (\\d+)-я (\\S+)ая [Дд]орожка ", " $1-ю $2ую дорожку "],
            ["^ [Дд]орожка ", " дорожку "],

            ["^ (\\S+)ая [Кк]оса ", " $1ую косу "],
            ["^ (\\S+)ая [Хх]орда ", " $1ую хорду "],

            ["^ [Дд]убл[её]р ", " дублёр "]
        ],
        "dative": [
            ["^ ([«\"])", " трасса $1"],

            ["^ (\\S+)ая [Аа]ллея ", " $1ой аллее "],
            ["^ (\\S+)ья [Аа]ллея ", " $1ьей аллее "],
            ["^ (\\S+)яя [Аа]ллея ", " $1ей аллее "],
            ["^ (\\d+)-я (\\S+)ая [Аа]ллея ", " $1-й $2ой аллее "],
            ["^ [Аа]ллея ", " аллее "],

            ["^ (\\S+)ая-(\\S+)ая [Уу]лица ", " $1ой-$2ой улице "],
            ["^ (\\S+)ая [Уу]лица ", " $1ой улице "],
            ["^ (\\S+)ья [Уу]лица ", " $1ьей улице "],
            ["^ (\\S+)яя [Уу]лица ", " $1ей улице "],
            ["^ (\\d+)-я [Уу]лица ", " $1-й улице "],
            ["^ (\\d+)-я (\\S+)ая [Уу]лица ", " $1-й $2ой улице "],
            ["^ (\\S+)ая (\\S+)ая [Уу]лица ", " $1ой $2ой улице "],
            ["^ (\\S+[вн])а [Уу]лица ", " $1ой улице "],
            ["^ (\\S+)ая (\\S+[вн])а [Уу]лица ", " $1ой $2ой улице "],
            ["^ Даньславля [Уу]лица ", " Даньславлей улице "],
            ["^ Добрыня [Уу]лица ", " Добрыней улице "],
            ["^ Людогоща [Уу]лица ", " Людогощей улице "],
            ["^ [Уу]лица ", " улице "],

            ["^ (\\d+)-я [Лл]иния ", " $1-й линии "],
            ["^ (\\d+)-(\\d+)-я [Лл]иния ", " $1-$2-й линии "],
            ["^ (\\S+)ая [Лл]иния ", " $1ой линии "],
            ["^ (\\S+)ья [Лл]иния ", " $1ьей линии "],
            ["^ (\\S+)яя [Лл]иния ", " $1ей линии "],
            ["^ (\\d+)-я (\\S+)ая [Лл]иния ", " $1-й $2ой линии "],
            ["^ [Лл]иния ", " линии "],

            ["^ (\\d+)-(\\d+)-я [Лл]инии ", " $1-$2-й линиям "],

            ["^ (\\S+)ая [Нн]абережная ", " $1ой набережной "],
            ["^ (\\S+)ья [Нн]абережная ", " $1ьей набережной "],
            ["^ (\\S+)яя [Нн]абережная ", " $1ей набережной "],
            ["^ (\\d+)-я (\\S+)ая [Нн]абережная ", " $1-й $2ой набережной "],
            ["^ [Нн]абережная ", " набережной "],

            ["^ (\\S+)ая [Пп]лощадь ", " $1ой площади "],
            ["^ (\\S+)ья [Пп]лощадь ", " $1ьей площади "],
            ["^ (\\S+)яя [Пп]лощадь ", " $1ей площади "],
            ["^ (\\S+[вн])а [Пп]лощадь ", " $1ой площади "],
            ["^ (\\d+)-я (\\S+)ая [Пп]лощадь ", " $1-й $2ой площади "],
            ["^ [Пп]лощадь ", " площади "],

            ["^ (\\S+)ая [Пп]росека ", " $1ой просеке "],
            ["^ (\\S+)ья [Пп]росека ", " $1ьей просеке "],
            ["^ (\\S+)яя [Пп]росека ", " $1ей просеке "],
            ["^ (\\d+)-я [Пп]росека ", " $1-й просеке "],
            ["^ [Пп]росека ", " просеке "],

            ["^ (\\S+)ая [Ээ]стакада ", " $1ой эстакаде "],
            ["^ (\\S+)ья [Ээ]стакада ", " $1ьей эстакаде "],
            ["^ (\\S+)яя [Ээ]стакада ", " $1ей эстакаде "],
            ["^ (\\d+)-я (\\S+)ая [Ээ]стакада ", " $1-й $2ой эстакаде "],
            ["^ [Ээ]стакада ", " эстакаде "],

            ["^ (\\S+)ая [Мм]агистраль ", " $1ой магистрали "],
            ["^ (\\S+)ья [Мм]агистраль ", " $1ьей магистрали "],
            ["^ (\\S+)яя [Мм]агистраль ", " $1ей магистрали "],
            ["^ (\\S+)ая (\\S+)ая [Мм]агистраль ", " $1ой $2ой магистрали "],
            ["^ (\\d+)-я (\\S+)ая [Мм]агистраль ", " $1-й $2ой магистрали "],
            ["^ [Мм]агистраль ", " магистрали "],

            ["^ (\\S+)ая [Рр]азвязка ", " $1ой развязке "],
            ["^ (\\S+)ья [Рр]азвязка ", " $1ьей развязке "],
            ["^ (\\S+)яя [Рр]азвязка ", " $1ей развязке "],
            ["^ (\\d+)-я (\\S+)ая [Рр]азвязка ", " $1-й $2ой развязке "],
            ["^ [Рр]азвязка ", " развязке "],

            ["^ (\\S+)ая [Тт]расса ", " $1ой трассе "],
            ["^ (\\S+)ья [Тт]расса ", " $1ьей трассе "],
            ["^ (\\S+)яя [Тт]расса ", " $1ей трассе "],
            ["^ (\\d+)-я (\\S+)ая [Тт]расса ", " $1-й $2ой трассе "],
            ["^ [Тт]расса ", " трассе "],

            ["^ (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2дороге "],
            ["^ (\\S+)ья ([Аа]вто)?[Дд]орога ", " $1ьей $2дороге "],
            ["^ (\\S+)яя ([Аа]вто)?[Дд]орога ", " $1ей $2дороге "],
            ["^ (\\S+)ая (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2ой $3дороге "],
            ["^ (\\d+)-я (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1-й $2ой $3дороге "],
            ["^ ([Аа]вто)?[Дд]орога ", " $1дороге "],

            ["^ (\\S+)ая [Дд]орожка ", " $1ой дорожке "],
            ["^ (\\S+)ья [Дд]орожка ", " $1ьей дорожке "],
            ["^ (\\S+)яя [Дд]орожка ", " $1ей дорожке "],
            ["^ (\\d+)-я (\\S+)ая [Дд]орожка ", " $1-й $2ой дорожке "],
            ["^ [Дд]орожка ", " дорожке "],

            ["^ (\\S+)во [Пп]оле ", " $1ву полю "],
            ["^ (\\S+)ая [Кк]оса ", " $1ой косе "],
            ["^ (\\S+)ая [Хх]орда ", " $1ой хорде "],
            ["^ (\\S+)[иоы]й [Пп]роток ", " $1ому протоку "],

            ["^ (\\S+н)ий [Бб]ульвар ", " $1ему бульвару "],
            ["^ (\\S+)[иоы]й [Бб]ульвар ", " $1ому бульвару "],
            ["^ (\\S+[иы]н) [Бб]ульвар ", " $1у бульвару "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Бб]ульвар ", " $1ому $2ему бульвару "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Бб]ульвар ", " $1ему $2ому бульвару "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Бб]ульвар ", " $1ому $2ому бульвару "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Бб]ульвар ", " $1ому $2у бульвару "],
            ["^ (\\d+)-й (\\S+н)ий [Бб]ульвар ", " $1-му $2ему бульвару "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Бб]ульвар ", " $1-му $2ому бульвару "],
            ["^ (\\d+)-й (\\S+[иы]н) [Бб]ульвар ", " $1-му $2у бульвару "],
            ["^ [Бб]ульвар ", " бульвару "],

            ["^ [Дд]убл[её]р ", " дублёру "],

            ["^ (\\S+н)ий [Зз]аезд ", " $1ему заезду "],
            ["^ (\\S+)[иоы]й [Зз]аезд ", " $1ому заезду "],
            ["^ (\\S+[еёо]в) [Зз]аезд ", " $1у заезду "],
            ["^ (\\S+[иы]н) [Зз]аезд ", " $1у заезду "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Зз]аезд ", " $1ому $2ему заезду "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Зз]аезд ", " $1ему $2ому заезду "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Зз]аезд ", " $1ому $2ому заезду "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Зз]аезд ", " $1ому $2у заезду "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Зз]аезд ", " $1ому $2у заезду "],
            ["^ (\\d+)-й (\\S+н)ий [Зз]аезд ", " $1-му $2ему заезду "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Зз]аезд ", " $1-му $2ому заезду "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Зз]аезд ", " $1-му $2у заезду "],
            ["^ (\\d+)-й (\\S+[иы]н) [Зз]аезд ", " $1-му $2у заезду "],
            ["^ [Зз]аезд ", " заезду "],

            ["^ (\\S+н)ий [Мм]ост ", " $1ему мосту "],
            ["^ (\\S+)[иоы]й [Мм]ост ", " $1ому мосту "],
            ["^ (\\S+[еёо]в) [Мм]ост ", " $1у мосту "],
            ["^ (\\S+[иы]н) [Мм]ост ", " $1у мосту "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Мм]ост ", " $1ому $2ему мосту "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Мм]ост ", " $1ему $2ому мосту "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Мм]ост ", " $1ому $2ому мосту "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Мм]ост ", " $1ому $2у мосту "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Мм]ост ", " $1ому $2у мосту "],
            ["^ (\\d+)-й [Мм]ост ", " $1-му мосту "],
            ["^ (\\d+)-й (\\S+н)ий [Мм]ост ", " $1-му $2ему мосту "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Мм]ост ", " $1-му $2ому мосту "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Мм]ост ", " $1-му $2у мосту "],
            ["^ (\\d+)-й (\\S+[иы]н) [Мм]ост ", " $1-му $2у мосту "],
            ["^ [Мм]ост ", " мосту "],

            ["^ (\\S+н)ий [Оо]бход ", " $1ему обходу "],
            ["^ (\\S+)[иоы]й [Оо]бход ", " $1ому обходу "],
            ["^ [Оо]бход ", " обходу "],

            ["^ (\\S+н)ий [Пп]арк ", " $1ему парку "],
            ["^ (\\S+)[иоы]й [Пп]арк ", " $1ому парку "],
            ["^ (\\S+[иы]н) [Пп]арк ", " $1у парку "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]арк ", " $1ому $2ему парку "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]арк ", " $1ему $2ому парку "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]арк ", " $1ому $2ому парку "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]арк ", " $1ому $2у парку "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]арк ", " $1-му $2ему парку "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]арк ", " $1-му $2ому парку "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]арк ", " $1-му $2у парку "],
            ["^ [Пп]арк ", " парку "],

            ["^ (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1ому-$2ому переулку "],
            ["^ (\\d+)-й (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1-му $2ому-$3ому переулку "],
            ["^ (\\S+н)ий [Пп]ереулок ", " $1ему переулку "],
            ["^ (\\S+)[иоы]й [Пп]ереулок ", " $1ому переулку "],
            ["^ (\\S+[еёо]в) [Пп]ереулок ", " $1у переулку "],
            ["^ (\\S+[иы]н) [Пп]ереулок ", " $1у переулку "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]ереулок ", " $1ому $2ему переулку "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]ереулок ", " $1ему $2ому переулку "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]ереулок ", " $1ому $2ому переулку "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]ереулок ", " $1ому $2у переулку "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]ереулок ", " $1ому $2у переулку "],
            ["^ (\\d+)-й [Пп]ереулок ", " $1-му переулку "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]ереулок ", " $1-му $2ему переулку "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]ереулок ", " $1-му $2ому переулку "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]ереулок ", " $1-му $2у переулку "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]ереулок ", " $1-му $2у переулку "],
            ["^ [Пп]ереулок ", " переулку "],

            ["^ [Пп]одъезд ", " подъезду "],

            ["^ (\\S+[еёо]в)-(\\S+)[иоы]й [Пп]роезд ", " $1у-$2ому проезду "],
            ["^ (\\S+н)ий [Пп]роезд ", " $1ему проезду "],
            ["^ (\\S+)[иоы]й [Пп]роезд ", " $1ому проезду "],
            ["^ (\\S+[еёо]в) [Пп]роезд ", " $1у проезду "],
            ["^ (\\S+[иы]н) [Пп]роезд ", " $1у проезду "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роезд ", " $1ому $2ему проезду "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1ему $2ому проезду "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1ому $2ому проезду "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]роезд ", " $1ому $2у проезду "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роезд ", " $1ому $2у проезду "],
            ["^ (\\d+)-й [Пп]роезд ", " $1-му проезду "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роезд ", " $1-му $2ему проезду "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роезд ", " $1-му $2ому проезду "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]роезд ", " $1-му $2у проезду "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роезд ", " $1-му $2у проезду "],
            ["^ (\\d+)-й (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1-му $2ему $3ому проезду "],
            ["^ (\\d+)-й (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1-му $2ому $3ому проезду "],
            ["^ [Пп]роезд ", " проезду "],

            ["^ (\\S+н)ий [Пп]роспект ", " $1ему проспекту "],
            ["^ (\\S+)[иоы]й [Пп]роспект ", " $1ому проспекту "],
            ["^ (\\S+[иы]н) [Пп]роспект ", " $1у проспекту "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роспект ", " $1ому $2ему проспекту "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роспект ", " $1ему $2ому проспекту "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роспект ", " $1ому $2ому проспекту "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роспект ", " $1ому $2у проспекту "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роспект ", " $1-му $2ему проспекту "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роспект ", " $1-му $2ому проспекту "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роспект ", " $1-му $2у проспекту "],
            ["^ [Пп]роспект ", " проспекту "],

            ["^ (\\S+н)ий [Пп]утепровод ", " $1ему путепроводу "],
            ["^ (\\S+)[иоы]й [Пп]утепровод ", " $1ому путепроводу "],
            ["^ (\\S+[иы]н) [Пп]утепровод ", " $1у путепроводу "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]утепровод ", " $1ому $2ему путепроводу "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]утепровод ", " $1ему $2ому путепроводу "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]утепровод ", " $1ому $2ому путепроводу "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]утепровод ", " $1ому $2у путепроводу "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]утепровод ", " $1-му $2ему путепроводу "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]утепровод ", " $1-му $2ому путепроводу "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]утепровод ", " $1-му $2у путепроводу "],
            ["^ [Пп]утепровод ", " путепроводу "],

            ["^ (\\S+н)ий [Сс]пуск ", " $1ему спуску "],
            ["^ (\\S+)[иоы]й [Сс]пуск ", " $1ому спуску "],
            ["^ (\\S+[еёо]в) [Сс]пуск ", " $1у спуску "],
            ["^ (\\S+[иы]н) [Сс]пуск ", " $1у спуску "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]пуск ", " $1ому $2ему спуску "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]пуск ", " $1ему $2ому спуску "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]пуск ", " $1ому $2ому спуску "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Сс]пуск ", " $1ому $2у спуску "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]пуск ", " $1ому $2у спуску "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]пуск ", " $1-му $2ему спуску "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]пуск ", " $1-му $2ому спуску "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Сс]пуск ", " $1-му $2у спуску "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]пуск ", " $1-му $2у спуску "],
            ["^ [Сс]пуск ", " спуску "],

            ["^ (\\S+н)ий [Сс]ъезд ", " $1ему съезду "],
            ["^ (\\S+)[иоы]й [Сс]ъезд ", " $1ому съезду "],
            ["^ (\\S+[иы]н) [Сс]ъезд ", " $1у съезду "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]ъезд ", " $1ому $2ему съезду "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]ъезд ", " $1ему $2ому съезду "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]ъезд ", " $1ому $2ому съезду "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]ъезд ", " $1ому $2у съезду "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]ъезд ", " $1-му $2ему съезду "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]ъезд ", " $1-му $2ому съезду "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]ъезд ", " $1-му $2у съезду "],
            ["^ [Сс]ъезд ", " съезду "],

            ["^ (\\S+н)ий [Тт][уо]ннель ", " $1ему тоннелю "],
            ["^ (\\S+)[иоы]й [Тт][уо]ннель ", " $1ому тоннелю "],
            ["^ (\\S+[иы]н) [Тт][уо]ннель ", " $1у тоннелю "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт][уо]ннель ", " $1ому $2ему тоннелю "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт][уо]ннель ", " $1ему $2ому тоннелю "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт][уо]ннель ", " $1ому $2ому тоннелю "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт][уо]ннель ", " $1ому $2у тоннелю "],
            ["^ (\\d+)-й (\\S+н)ий [Тт][уо]ннель ", " $1-му $2ему тоннелю "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт][уо]ннель ", " $1-му $2ому тоннелю "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт][уо]ннель ", " $1-му $2у тоннелю "],
            ["^ [Тт][уо]ннель ", " тоннелю "],

            ["^ (\\S+н)ий [Тт]ракт ", " $1ему тракту "],
            ["^ (\\S+)[иоы]й [Тт]ракт ", " $1ому тракту "],
            ["^ (\\S+[еёо]в) [Тт]ракт ", " $1у тракту "],
            ["^ (\\S+[иы]н) [Тт]ракт ", " $1у тракту "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]ракт ", " $1ому $2ему тракту "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]ракт ", " $1ему $2ому тракту "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]ракт ", " $1ому $2ому тракту "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]ракт ", " $1ому $2у тракту "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]ракт ", " $1ому $2у тракту "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]ракт ", " $1-му $2ему тракту "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]ракт ", " $1-му $2ому тракту "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]ракт ", " $1-му $2у тракту "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]ракт ", " $1-му $2у тракту "],
            ["^ [Тт]ракт ", " тракту "],

            ["^ (\\S+н)ий [Тт]упик ", " $1ему тупику "],
            ["^ (\\S+)[иоы]й [Тт]упик ", " $1ому тупику "],
            ["^ (\\S+[еёо]в) [Тт]упик ", " $1у тупику "],
            ["^ (\\S+[иы]н) [Тт]упик ", " $1у тупику "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]упик ", " $1ому $2ему тупику "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]упик ", " $1ему $2ому тупику "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]упик ", " $1ому $2ому тупику "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]упик ", " $1ому $2у тупику "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]упик ", " $1ому $2у тупику "],
            ["^ (\\d+)-й [Тт]упик ", " $1-му тупику "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]упик ", " $1-му $2ему тупику "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]упик ", " $1-му $2ому тупику "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]упик ", " $1-му $2у тупику "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]упик ", " $1-му $2у тупику "],
            ["^ [Тт]упик ", " тупику "],

            ["^ (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1му $2кольцу "],
            ["^ (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1му $2кольцу "],
            ["^ (\\S+[ео])е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1му $2му $3кольцу "],
            ["^ (\\S+ье) (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1му $2му $3кольцу "],
            ["^ (\\d+)-е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1-му $2му $3кольцу "],
            ["^ (\\d+)-е (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1-му $2му $3кольцу "],
            ["^ ([Пп]олу)?[Кк]ольцо ", " $1кольцу "],

            ["^ (\\S+[ео])е [Шш]оссе ", " $1му шоссе "],
            ["^ (\\S+ье) [Шш]оссе ", " $1му шоссе "],
            ["^ (\\S+[ео])е (\\S+[ео])е [Шш]оссе ", " $1му $2му шоссе "],
            ["^ (\\S+ье) (\\S+[ео])е [Шш]оссе ", " $1му $2му шоссе "],
            ["^ (\\d+)-е (\\S+[ео])е [Шш]оссе ", " $1-му $2му шоссе "],
            ["^ (\\d+)-е (\\S+ье) [Шш]оссе ", " $1-му $2му шоссе "],

            [" ([Тт])ретому ", " $1ретьему "],
            ["([жч])ому ", "$1ьему "],
            ["([жч])ой ", "$1ей "]
        ],
        "genitive": [
            ["^ ([«\"])", " трасса $1"],

            ["^ (\\S+)ая [Аа]ллея ", " $1ой аллеи "],
            ["^ (\\S+)ья [Аа]ллея ", " $1ьей аллеи "],
            ["^ (\\S+)яя [Аа]ллея ", " $1ей аллеи "],
            ["^ (\\d+)-я (\\S+)ая [Аа]ллея ", " $1-й $2ой аллеи "],
            ["^ [Аа]ллея ", " аллеи "],

            ["^ (\\S+)ая-(\\S+)ая [Уу]лица ", " $1ой-$2ой улицы "],
            ["^ (\\S+)ая [Уу]лица ", " $1ой улицы "],
            ["^ (\\S+)ья [Уу]лица ", " $1ьей улицы "],
            ["^ (\\S+)яя [Уу]лица ", " $1ей улицы "],
            ["^ (\\d+)-я [Уу]лица ", " $1-й улицы "],
            ["^ (\\d+)-я (\\S+)ая [Уу]лица ", " $1-й $2ой улицы "],
            ["^ (\\S+)ая (\\S+)ая [Уу]лица ", " $1ой $2ой улицы "],
            ["^ (\\S+[вн])а [Уу]лица ", " $1ой улицы "],
            ["^ (\\S+)ая (\\S+[вн])а [Уу]лица ", " $1ой $2ой улицы "],
            ["^ Даньславля [Уу]лица ", " Даньславлей улицы "],
            ["^ Добрыня [Уу]лица ", " Добрыней улицы "],
            ["^ Людогоща [Уу]лица ", " Людогощей улицы "],
            ["^ [Уу]лица ", " улицы "],

            ["^ (\\d+)-я [Лл]иния ", " $1-й линии "],
            ["^ (\\d+)-(\\d+)-я [Лл]иния ", " $1-$2-й линии "],
            ["^ (\\S+)ая [Лл]иния ", " $1ой линии "],
            ["^ (\\S+)ья [Лл]иния ", " $1ьей линии "],
            ["^ (\\S+)яя [Лл]иния ", " $1ей линии "],
            ["^ (\\d+)-я (\\S+)ая [Лл]иния ", " $1-й $2ой линии "],
            ["^ [Лл]иния ", " линии "],

            ["^ (\\d+)-(\\d+)-я [Лл]инии ", " $1-$2-й линий "],

            ["^ (\\S+)ая [Нн]абережная ", " $1ой набережной "],
            ["^ (\\S+)ья [Нн]абережная ", " $1ьей набережной "],
            ["^ (\\S+)яя [Нн]абережная ", " $1ей набережной "],
            ["^ (\\d+)-я (\\S+)ая [Нн]абережная ", " $1-й $2ой набережной "],
            ["^ [Нн]абережная ", " набережной "],

            ["^ (\\S+)ая [Пп]лощадь ", " $1ой площади "],
            ["^ (\\S+)ья [Пп]лощадь ", " $1ьей площади "],
            ["^ (\\S+)яя [Пп]лощадь ", " $1ей площади "],
            ["^ (\\S+[вн])а [Пп]лощадь ", " $1ой площади "],
            ["^ (\\d+)-я (\\S+)ая [Пп]лощадь ", " $1-й $2ой площади "],
            ["^ [Пп]лощадь ", " площади "],

            ["^ (\\S+)ая [Пп]росека ", " $1ой просеки "],
            ["^ (\\S+)ья [Пп]росека ", " $1ьей просеки "],
            ["^ (\\S+)яя [Пп]росека ", " $1ей просеки "],
            ["^ (\\d+)-я [Пп]росека ", " $1-й просеки "],
            ["^ [Пп]росека ", " просеки "],

            ["^ (\\S+)ая [Ээ]стакада ", " $1ой эстакады "],
            ["^ (\\S+)ья [Ээ]стакада ", " $1ьей эстакады "],
            ["^ (\\S+)яя [Ээ]стакада ", " $1ей эстакады "],
            ["^ (\\d+)-я (\\S+)ая [Ээ]стакада ", " $1-й $2ой эстакады "],
            ["^ [Ээ]стакада ", " эстакады "],

            ["^ (\\S+)ая [Мм]агистраль ", " $1ой магистрали "],
            ["^ (\\S+)ья [Мм]агистраль ", " $1ьей магистрали "],
            ["^ (\\S+)яя [Мм]агистраль ", " $1ей магистрали "],
            ["^ (\\S+)ая (\\S+)ая [Мм]агистраль ", " $1ой $2ой магистрали "],
            ["^ (\\d+)-я (\\S+)ая [Мм]агистраль ", " $1-й $2ой магистрали "],
            ["^ [Мм]агистраль ", " магистрали "],

            ["^ (\\S+)ая [Рр]азвязка ", " $1ой развязки "],
            ["^ (\\S+)ья [Рр]азвязка ", " $1ьей развязки "],
            ["^ (\\S+)яя [Рр]азвязка ", " $1ей развязки "],
            ["^ (\\d+)-я (\\S+)ая [Рр]азвязка ", " $1-й $2ой развязки "],
            ["^ [Рр]азвязка ", " развязки "],

            ["^ (\\S+)ая [Тт]расса ", " $1ой трассы "],
            ["^ (\\S+)ья [Тт]расса ", " $1ьей трассы "],
            ["^ (\\S+)яя [Тт]расса ", " $1ей трассы "],
            ["^ (\\d+)-я (\\S+)ая [Тт]расса ", " $1-й $2ой трассы "],
            ["^ [Тт]расса ", " трассы "],

            ["^ (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2дороги "],
            ["^ (\\S+)ья ([Аа]вто)?[Дд]орога ", " $1ьей $2дороги "],
            ["^ (\\S+)яя ([Аа]вто)?[Дд]орога ", " $1ей $2дороги "],
            ["^ (\\S+)ая (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2ой $3дороги "],
            ["^ (\\d+)-я (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1-й $2ой $3дороги "],
            ["^ ([Аа]вто)?[Дд]орога ", " $1дороги "],

            ["^ (\\S+)ая [Дд]орожка ", " $1ой дорожки "],
            ["^ (\\S+)ья [Дд]орожка ", " $1ьей дорожки "],
            ["^ (\\S+)яя [Дд]орожка ", " $1ей дорожки "],
            ["^ (\\d+)-я (\\S+)ая [Дд]орожка ", " $1-й $2ой дорожки "],
            ["^ [Дд]орожка ", " дорожки "],

            ["^ (\\S+)во [Пп]оле ", " $1ва поля "],
            ["^ (\\S+)ая [Кк]оса ", " $1ой косы "],
            ["^ (\\S+)ая [Хх]орда ", " $1ой хорды "],
            ["^ (\\S+)[иоы]й [Пп]роток ", " $1ого протока "],

            ["^ (\\S+н)ий [Бб]ульвар ", " $1его бульвара "],
            ["^ (\\S+)[иоы]й [Бб]ульвар ", " $1ого бульвара "],
            ["^ (\\S+[иы]н) [Бб]ульвар ", " $1ого бульвара "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Бб]ульвар ", " $1ого $2его бульвара "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Бб]ульвар ", " $1его $2ого бульвара "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Бб]ульвар ", " $1ого $2ого бульвара "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Бб]ульвар ", " $1ого $2ого бульвара "],
            ["^ (\\d+)-й (\\S+н)ий [Бб]ульвар ", " $1-го $2его бульвара "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Бб]ульвар ", " $1-го $2ого бульвара "],
            ["^ (\\d+)-й (\\S+[иы]н) [Бб]ульвар ", " $1-го $2ого бульвара "],
            ["^ [Бб]ульвар ", " бульвара "],

            ["^ [Дд]убл[её]р ", " дублёра "],

            ["^ (\\S+н)ий [Зз]аезд ", " $1его заезда "],
            ["^ (\\S+)[иоы]й [Зз]аезд ", " $1ого заезда "],
            ["^ (\\S+[еёо]в) [Зз]аезд ", " $1а заезда "],
            ["^ (\\S+[иы]н) [Зз]аезд ", " $1а заезда "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Зз]аезд ", " $1ого $2его заезда "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Зз]аезд ", " $1его $2ого заезда "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Зз]аезд ", " $1ого $2ого заезда "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Зз]аезд ", " $1ого $2а заезда "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Зз]аезд ", " $1ого $2а заезда "],
            ["^ (\\d+)-й (\\S+н)ий [Зз]аезд ", " $1-го $2его заезда "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Зз]аезд ", " $1-го $2ого заезда "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Зз]аезд ", " $1-го $2а заезда "],
            ["^ (\\d+)-й (\\S+[иы]н) [Зз]аезд ", " $1-го $2а заезда "],
            ["^ [Зз]аезд ", " заезда "],

            ["^ (\\S+н)ий [Мм]ост ", " $1его моста "],
            ["^ (\\S+)[иоы]й [Мм]ост ", " $1ого моста "],
            ["^ (\\S+[еёо]в) [Мм]ост ", " $1а моста "],
            ["^ (\\S+[иы]н) [Мм]ост ", " $1а моста "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Мм]ост ", " $1ого $2его моста "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Мм]ост ", " $1его $2ого моста "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Мм]ост ", " $1ого $2ого моста "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Мм]ост ", " $1ого $2а моста "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Мм]ост ", " $1ого $2а моста "],
            ["^ (\\d+)-й [Мм]ост ", " $1-го моста "],
            ["^ (\\d+)-й (\\S+н)ий [Мм]ост ", " $1-го $2его моста "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Мм]ост ", " $1-го $2ого моста "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Мм]ост ", " $1-го $2а моста "],
            ["^ (\\d+)-й (\\S+[иы]н) [Мм]ост ", " $1-го $2а моста "],
            ["^ [Мм]ост ", " моста "],

            ["^ (\\S+н)ий [Оо]бход ", " $1его обхода "],
            ["^ (\\S+)[иоы]й [Оо]бход ", " $1ого обхода "],
            ["^ [Оо]бход ", " обхода "],

            ["^ (\\S+н)ий [Пп]арк ", " $1его парка "],
            ["^ (\\S+)[иоы]й [Пп]арк ", " $1ого парка "],
            ["^ (\\S+[иы]н) [Пп]арк ", " $1ого парка "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]арк ", " $1ого $2его парка "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]арк ", " $1его $2ого парка "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]арк ", " $1ого $2ого парка "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]арк ", " $1ого $2ого парка "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]арк ", " $1-го $2его парка "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]арк ", " $1-го $2ого парка "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]арк ", " $1-го $2ого парка "],
            ["^ [Пп]арк ", " парка "],

            ["^ (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1ого-$2ого переулка "],
            ["^ (\\d+)-й (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1-го $2ого-$3ого переулка "],
            ["^ (\\S+н)ий [Пп]ереулок ", " $1его переулка "],
            ["^ (\\S+)[иоы]й [Пп]ереулок ", " $1ого переулка "],
            ["^ (\\S+[еёо]в) [Пп]ереулок ", " $1а переулка "],
            ["^ (\\S+[иы]н) [Пп]ереулок ", " $1а переулка "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]ереулок ", " $1ого $2его переулка "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]ереулок ", " $1его $2ого переулка "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]ереулок ", " $1ого $2ого переулка "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]ереулок ", " $1ого $2а переулка "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]ереулок ", " $1ого $2а переулка "],
            ["^ (\\d+)-й [Пп]ереулок ", " $1-го переулка "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]ереулок ", " $1-го $2его переулка "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]ереулок ", " $1-го $2ого переулка "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]ереулок ", " $1-го $2а переулка "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]ереулок ", " $1-го $2а переулка "],
            ["^ [Пп]ереулок ", " переулка "],

            ["^ [Пп]одъезд ", " подъезда "],

            ["^ (\\S+[еёо]в)-(\\S+)[иоы]й [Пп]роезд ", " $1а-$2ого проезда "],
            ["^ (\\S+н)ий [Пп]роезд ", " $1его проезда "],
            ["^ (\\S+)[иоы]й [Пп]роезд ", " $1ого проезда "],
            ["^ (\\S+[еёо]в) [Пп]роезд ", " $1а проезда "],
            ["^ (\\S+[иы]н) [Пп]роезд ", " $1а проезда "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роезд ", " $1ого $2его проезда "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1его $2ого проезда "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1ого $2ого проезда "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]роезд ", " $1ого $2а проезда "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роезд ", " $1ого $2а проезда "],
            ["^ (\\d+)-й [Пп]роезд ", " $1-го проезда "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роезд ", " $1-го $2его проезда "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роезд ", " $1-го $2ого проезда "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]роезд ", " $1-го $2а проезда "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роезд ", " $1-го $2а проезда "],
            ["^ (\\d+)-й (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1-го $2его $3ого проезда "],
            ["^ (\\d+)-й (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1-го $2ого $3ого проезда "],
            ["^ [Пп]роезд ", " проезда "],

            ["^ (\\S+н)ий [Пп]роспект ", " $1его проспекта "],
            ["^ (\\S+)[иоы]й [Пп]роспект ", " $1ого проспекта "],
            ["^ (\\S+[иы]н) [Пп]роспект ", " $1ого проспекта "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роспект ", " $1ого $2его проспекта "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роспект ", " $1его $2ого проспекта "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роспект ", " $1ого $2ого проспекта "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роспект ", " $1ого $2ого проспекта "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роспект ", " $1-го $2его проспекта "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роспект ", " $1-го $2ого проспекта "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роспект ", " $1-го $2ого проспекта "],
            ["^ [Пп]роспект ", " проспекта "],

            ["^ (\\S+н)ий [Пп]утепровод ", " $1его путепровода "],
            ["^ (\\S+)[иоы]й [Пп]утепровод ", " $1ого путепровода "],
            ["^ (\\S+[иы]н) [Пп]утепровод ", " $1ого путепровода "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]утепровод ", " $1ого $2его путепровода "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]утепровод ", " $1его $2ого путепровода "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]утепровод ", " $1ого $2ого путепровода "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]утепровод ", " $1ого $2ого путепровода "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]утепровод ", " $1-го $2его путепровода "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]утепровод ", " $1-го $2ого путепровода "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]утепровод ", " $1-го $2ого путепровода "],
            ["^ [Пп]утепровод ", " путепровода "],

            ["^ (\\S+н)ий [Сс]пуск ", " $1его спуска "],
            ["^ (\\S+)[иоы]й [Сс]пуск ", " $1ого спуска "],
            ["^ (\\S+[еёо]в) [Сс]пуск ", " $1а спуска "],
            ["^ (\\S+[иы]н) [Сс]пуск ", " $1а спуска "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]пуск ", " $1ого $2его спуска "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]пуск ", " $1его $2ого спуска "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]пуск ", " $1ого $2ого спуска "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Сс]пуск ", " $1ого $2а спуска "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]пуск ", " $1ого $2а спуска "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]пуск ", " $1-го $2его спуска "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]пуск ", " $1-го $2ого спуска "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Сс]пуск ", " $1-го $2а спуска "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]пуск ", " $1-го $2а спуска "],
            ["^ [Сс]пуск ", " спуска "],

            ["^ (\\S+н)ий [Сс]ъезд ", " $1его съезда "],
            ["^ (\\S+)[иоы]й [Сс]ъезд ", " $1ого съезда "],
            ["^ (\\S+[иы]н) [Сс]ъезд ", " $1ого съезда "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]ъезд ", " $1ого $2его съезда "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]ъезд ", " $1его $2ого съезда "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]ъезд ", " $1ого $2ого съезда "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]ъезд ", " $1ого $2ого съезда "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]ъезд ", " $1-го $2его съезда "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]ъезд ", " $1-го $2ого съезда "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]ъезд ", " $1-го $2ого съезда "],
            ["^ [Сс]ъезд ", " съезда "],

            ["^ (\\S+н)ий [Тт][уо]ннель ", " $1его тоннеля "],
            ["^ (\\S+)[иоы]й [Тт][уо]ннель ", " $1ого тоннеля "],
            ["^ (\\S+[иы]н) [Тт][уо]ннель ", " $1ого тоннеля "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт][уо]ннель ", " $1ого $2его тоннеля "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт][уо]ннель ", " $1его $2ого тоннеля "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт][уо]ннель ", " $1ого $2ого тоннеля "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт][уо]ннель ", " $1ого $2ого тоннеля "],
            ["^ (\\d+)-й (\\S+н)ий [Тт][уо]ннель ", " $1-го $2его тоннеля "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт][уо]ннель ", " $1-го $2ого тоннеля "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт][уо]ннель ", " $1-го $2ого тоннеля "],
            ["^ [Тт][уо]ннель ", " тоннеля "],

            ["^ (\\S+н)ий [Тт]ракт ", " $1ем тракта "],
            ["^ (\\S+)[иоы]й [Тт]ракт ", " $1ого тракта "],
            ["^ (\\S+[еёо]в) [Тт]ракт ", " $1а тракта "],
            ["^ (\\S+[иы]н) [Тт]ракт ", " $1а тракта "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]ракт ", " $1ого $2его тракта "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]ракт ", " $1его $2ого тракта "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]ракт ", " $1ого $2ого тракта "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]ракт ", " $1ого $2а тракта "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]ракт ", " $1ого $2а тракта "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]ракт ", " $1-го $2его тракта "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]ракт ", " $1-го $2ого тракта "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]ракт ", " $1-го $2а тракта "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]ракт ", " $1-го $2а тракта "],
            ["^ [Тт]ракт ", " тракта "],

            ["^ (\\S+н)ий [Тт]упик ", " $1его тупика "],
            ["^ (\\S+)[иоы]й [Тт]упик ", " $1ого тупика "],
            ["^ (\\S+[еёо]в) [Тт]упик ", " $1а тупика "],
            ["^ (\\S+[иы]н) [Тт]упик ", " $1а тупика "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]упик ", " $1ого $2его тупика "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]упик ", " $1его $2ого тупика "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]упик ", " $1ого $2ого тупика "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]упик ", " $1ого $2а тупика "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]упик ", " $1ого $2а тупика "],
            ["^ (\\d+)-й [Тт]упик ", " $1-го тупика "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]упик ", " $1-го $2его тупика "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]упик ", " $1-го $2ого тупика "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]упик ", " $1-го $2а тупика "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]упик ", " $1-го $2а тупика "],
            ["^ [Тт]упик ", " тупика "],

            ["^ (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1го $2кольца "],
            ["^ (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1го $2кольца "],
            ["^ (\\S+[ео])е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1го $2го $3кольца "],
            ["^ (\\S+ье) (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1го $2го $3кольца "],
            ["^ (\\d+)-е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1-го $2го $3кольца "],
            ["^ (\\d+)-е (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1-го $2го $3кольца "],
            ["^ ([Пп]олу)?[Кк]ольцо ", " $1кольца "],

            ["^ (\\S+[ео])е [Шш]оссе ", " $1го шоссе "],
            ["^ (\\S+ье) [Шш]оссе ", " $1го шоссе "],
            ["^ (\\S+[ео])е (\\S+[ео])е [Шш]оссе ", " $1го $2го шоссе "],
            ["^ (\\S+ье) (\\S+[ео])е [Шш]оссе ", " $1го $2го шоссе "],
            ["^ (\\d+)-е (\\S+[ео])е [Шш]оссе ", " $1-го $2го шоссе "],
            ["^ (\\d+)-е (\\S+ье) [Шш]оссе ", " $1-го $2го шоссе "],

            [" ([Тт])ретого ", " $1ретьего "],
            ["([жч])ого ", "$1ьего "]
        ],
        "prepositional": [
            ["^ ([«\"])", " трасса $1"],

            ["^ (\\S+)ая [Аа]ллея ", " $1ой аллее "],
            ["^ (\\S+)ья [Аа]ллея ", " $1ьей аллее "],
            ["^ (\\S+)яя [Аа]ллея ", " $1ей аллее "],
            ["^ (\\d+)-я (\\S+)ая [Аа]ллея ", " $1-й $2ой аллее "],
            ["^ [Аа]ллея ", " аллее "],

            ["^ (\\S+)ая-(\\S+)ая [Уу]лица ", " $1ой-$2ой улице "],
            ["^ (\\S+)ая [Уу]лица ", " $1ой улице "],
            ["^ (\\S+)ья [Уу]лица ", " $1ьей улице "],
            ["^ (\\S+)яя [Уу]лица ", " $1ей улице "],
            ["^ (\\d+)-я [Уу]лица ", " $1-й улице "],
            ["^ (\\d+)-я (\\S+)ая [Уу]лица ", " $1-й $2ой улице "],
            ["^ (\\S+)ая (\\S+)ая [Уу]лица ", " $1ой $2ой улице "],
            ["^ (\\S+[вн])а [Уу]лица ", " $1ой улице "],
            ["^ (\\S+)ая (\\S+[вн])а [Уу]лица ", " $1ой $2ой улице "],
            ["^ Даньславля [Уу]лица ", " Даньславлей улице "],
            ["^ Добрыня [Уу]лица ", " Добрыней улице "],
            ["^ Людогоща [Уу]лица ", " Людогощей улице "],
            ["^ [Уу]лица ", " улице "],

            ["^ (\\d+)-я [Лл]иния ", " $1-й линии "],
            ["^ (\\d+)-(\\d+)-я [Лл]иния ", " $1-$2-й линии "],
            ["^ (\\S+)ая [Лл]иния ", " $1ой линии "],
            ["^ (\\S+)ья [Лл]иния ", " $1ьей линии "],
            ["^ (\\S+)яя [Лл]иния ", " $1ей линии "],
            ["^ (\\d+)-я (\\S+)ая [Лл]иния ", " $1-й $2ой линии "],
            ["^ [Лл]иния ", " линии "],

            ["^ (\\d+)-(\\d+)-я [Лл]инии ", " $1-$2-й линиях "],

            ["^ (\\S+)ая [Нн]абережная ", " $1ой набережной "],
            ["^ (\\S+)ья [Нн]абережная ", " $1ьей набережной "],
            ["^ (\\S+)яя [Нн]абережная ", " $1ей набережной "],
            ["^ (\\d+)-я (\\S+)ая [Нн]абережная ", " $1-й $2ой набережной "],
            ["^ [Нн]абережная ", " набережной "],

            ["^ (\\S+)ая [Пп]лощадь ", " $1ой площади "],
            ["^ (\\S+)ья [Пп]лощадь ", " $1ьей площади "],
            ["^ (\\S+)яя [Пп]лощадь ", " $1ей площади "],
            ["^ (\\S+[вн])а [Пп]лощадь ", " $1ой площади "],
            ["^ (\\d+)-я (\\S+)ая [Пп]лощадь ", " $1-й $2ой площади "],
            ["^ [Пп]лощадь ", " площади "],

            ["^ (\\S+)ая [Пп]росека ", " $1ой просеке "],
            ["^ (\\S+)ья [Пп]росека ", " $1ьей просеке "],
            ["^ (\\S+)яя [Пп]росека ", " $1ей просеке "],
            ["^ (\\d+)-я [Пп]росека ", " $1-й просеке "],
            ["^ [Пп]росека ", " просеке "],

            ["^ (\\S+)ая [Ээ]стакада ", " $1ой эстакаде "],
            ["^ (\\S+)ья [Ээ]стакада ", " $1ьей эстакаде "],
            ["^ (\\S+)яя [Ээ]стакада ", " $1ей эстакаде "],
            ["^ (\\d+)-я (\\S+)ая [Ээ]стакада ", " $1-й $2ой эстакаде "],
            ["^ [Ээ]стакада ", " эстакаде "],

            ["^ (\\S+)ая [Мм]агистраль ", " $1ой магистрали "],
            ["^ (\\S+)ья [Мм]агистраль ", " $1ьей магистрали "],
            ["^ (\\S+)яя [Мм]агистраль ", " $1ей магистрали "],
            ["^ (\\S+)ая (\\S+)ая [Мм]агистраль ", " $1ой $2ой магистрали "],
            ["^ (\\d+)-я (\\S+)ая [Мм]агистраль ", " $1-й $2ой магистрали "],
            ["^ [Мм]агистраль ", " магистрали "],

            ["^ (\\S+)ая [Рр]азвязка ", " $1ой развязке "],
            ["^ (\\S+)ья [Рр]азвязка ", " $1ьей развязке "],
            ["^ (\\S+)яя [Рр]азвязка ", " $1ей развязке "],
            ["^ (\\d+)-я (\\S+)ая [Рр]азвязка ", " $1-й $2ой развязке "],
            ["^ [Рр]азвязка ", " развязке "],

            ["^ (\\S+)ая [Тт]расса ", " $1ой трассе "],
            ["^ (\\S+)ья [Тт]расса ", " $1ьей трассе "],
            ["^ (\\S+)яя [Тт]расса ", " $1ей трассе "],
            ["^ (\\d+)-я (\\S+)ая [Тт]расса ", " $1-й $2ой трассе "],
            ["^ [Тт]расса ", " трассе "],

            ["^ (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2дороге "],
            ["^ (\\S+)ья ([Аа]вто)?[Дд]орога ", " $1ьей $2дороге "],
            ["^ (\\S+)яя ([Аа]вто)?[Дд]орога ", " $1ей $2дороге "],
            ["^ (\\S+)ая (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1ой $2ой $3дороге "],
            ["^ (\\d+)-я (\\S+)ая ([Аа]вто)?[Дд]орога ", " $1-й $2ой $3дороге "],
            ["^ ([Аа]вто)?[Дд]орога ", " $1дороге "],

            ["^ (\\S+)ая [Дд]орожка ", " $1ой дорожке "],
            ["^ (\\S+)ья [Дд]орожка ", " $1ьей дорожке "],
            ["^ (\\S+)яя [Дд]орожка ", " $1ей дорожке "],
            ["^ (\\d+)-я (\\S+)ая [Дд]орожка ", " $1-й $2ой дорожке "],
            ["^ [Дд]орожка ", " дорожке "],

            ["^ (\\S+)во [Пп]оле ", " $1вом поле "],
            ["^ (\\S+)ая [Кк]оса ", " $1ой косе "],
            ["^ (\\S+)ая [Хх]орда ", " $1ой хорде "],
            ["^ (\\S+)[иоы]й [Пп]роток ", " $1ом протоке "],

            ["^ (\\S+н)ий [Бб]ульвар ", " $1ем бульваре "],
            ["^ (\\S+)[иоы]й [Бб]ульвар ", " $1ом бульваре "],
            ["^ (\\S+[иы]н) [Бб]ульвар ", " $1ом бульваре "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Бб]ульвар ", " $1ом $2ем бульваре "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Бб]ульвар ", " $1ем $2ом бульваре "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Бб]ульвар ", " $1ом $2ом бульваре "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Бб]ульвар ", " $1ом $2ом бульваре "],
            ["^ (\\d+)-й (\\S+н)ий [Бб]ульвар ", " $1-м $2ем бульваре "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Бб]ульвар ", " $1-м $2ом бульваре "],
            ["^ (\\d+)-й (\\S+[иы]н) [Бб]ульвар ", " $1-м $2ом бульваре "],
            ["^ [Бб]ульвар ", " бульваре "],

            ["^ [Дд]убл[её]р ", " дублёре "],

            ["^ (\\S+н)ий [Зз]аезд ", " $1ем заезде "],
            ["^ (\\S+)[иоы]й [Зз]аезд ", " $1ом заезде "],
            ["^ (\\S+[еёо]в) [Зз]аезд ", " $1ом заезде "],
            ["^ (\\S+[иы]н) [Зз]аезд ", " $1ом заезде "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Зз]аезд ", " $1ом $2ем заезде "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Зз]аезд ", " $1ем $2ом заезде "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Зз]аезд ", " $1ом $2ом заезде "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Зз]аезд ", " $1ом $2ом заезде "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Зз]аезд ", " $1ом $2ом заезде "],
            ["^ (\\d+)-й (\\S+н)ий [Зз]аезд ", " $1-м $2ем заезде "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Зз]аезд ", " $1-м $2ом заезде "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Зз]аезд ", " $1-м $2ом заезде "],
            ["^ (\\d+)-й (\\S+[иы]н) [Зз]аезд ", " $1-м $2ом заезде "],
            ["^ [Зз]аезд ", " заезде "],

            ["^ (\\S+н)ий [Мм]ост ", " $1ем мосту "],
            ["^ (\\S+)[иоы]й [Мм]ост ", " $1ом мосту "],
            ["^ (\\S+[еёо]в) [Мм]ост ", " $1ом мосту "],
            ["^ (\\S+[иы]н) [Мм]ост ", " $1ом мосту "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Мм]ост ", " $1ом $2ем мосту "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Мм]ост ", " $1ем $2ом мосту "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Мм]ост ", " $1ом $2ом мосту "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Мм]ост ", " $1ом $2ом мосту "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Мм]ост ", " $1ом $2ом мосту "],
            ["^ (\\d+)-й [Мм]ост ", " $1-м мосту "],
            ["^ (\\d+)-й (\\S+н)ий [Мм]ост ", " $1-м $2ем мосту "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Мм]ост ", " $1-м $2ом мосту "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Мм]ост ", " $1-м $2ом мосту "],
            ["^ (\\d+)-й (\\S+[иы]н) [Мм]ост ", " $1-м $2ом мосту "],
            ["^ [Мм]ост ", " мосту "],

            ["^ (\\S+н)ий [Оо]бход ", " $1ем обходе "],
            ["^ (\\S+)[иоы]й [Оо]бход ", " $1ом обходе "],
            ["^ [Оо]бход ", " обходе "],

            ["^ (\\S+н)ий [Пп]арк ", " $1ем парке "],
            ["^ (\\S+)[иоы]й [Пп]арк ", " $1ом парке "],
            ["^ (\\S+[иы]н) [Пп]арк ", " $1ом парке "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]арк ", " $1ом $2ем парке "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]арк ", " $1ем $2ом парке "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]арк ", " $1ом $2ом парке "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]арк ", " $1ом $2ом парке "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]арк ", " $1-м $2ем парке "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]арк ", " $1-м $2ом парке "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]арк ", " $1-м $2ом парке "],
            ["^ [Пп]арк ", " парке "],

            ["^ (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1ом-$2ом переулке "],
            ["^ (\\d+)-й (\\S+)[иоы]й-(\\S+)[иоы]й [Пп]ереулок ", " $1-м $2ом-$3ом переулке "],
            ["^ (\\S+н)ий [Пп]ереулок ", " $1ем переулке "],
            ["^ (\\S+)[иоы]й [Пп]ереулок ", " $1ом переулке "],
            ["^ (\\S+[еёо]в) [Пп]ереулок ", " $1ом переулке "],
            ["^ (\\S+[иы]н) [Пп]ереулок ", " $1ом переулке "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]ереулок ", " $1ом $2ем переулке "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]ереулок ", " $1ем $2ом переулке "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]ереулок ", " $1ом $2ом переулке "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]ереулок ", " $1ом $2ом переулке "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]ереулок ", " $1ом $2ом переулке "],
            ["^ (\\d+)-й [Пп]ереулок ", " $1-м переулке "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]ереулок ", " $1-м $2ем переулке "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]ереулок ", " $1-м $2ом переулке "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]ереулок ", " $1-м $2ом переулке "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]ереулок ", " $1-м $2ом переулке "],
            ["^ [Пп]ереулок ", " переулке "],

            ["^ [Пп]одъезд ", " подъезде "],

            ["^ (\\S+[еёо]в)-(\\S+)[иоы]й [Пп]роезд ", " $1ом-$2ом проезде "],
            ["^ (\\S+н)ий [Пп]роезд ", " $1ем проезде "],
            ["^ (\\S+)[иоы]й [Пп]роезд ", " $1ом проезде "],
            ["^ (\\S+[еёо]в) [Пп]роезд ", " $1ом проезде "],
            ["^ (\\S+[иы]н) [Пп]роезд ", " $1ом проезде "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роезд ", " $1ом $2ем проезде "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1ем $2ом проезде "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1ом $2ом проезде "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Пп]роезд ", " $1ом $2ом проезде "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роезд ", " $1ом $2ом проезде "],
            ["^ (\\d+)-й [Пп]роезд ", " $1-м проезде "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роезд ", " $1-м $2ем проезде "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роезд ", " $1-м $2ом проезде "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Пп]роезд ", " $1-м $2ом проезде "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роезд ", " $1-м $2ом проезде "],
            ["^ (\\d+)-й (\\S+н)ий (\\S+)[иоы]й [Пп]роезд ", " $1-м $2ем $3ом проезде "],
            ["^ (\\d+)-й (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роезд ", " $1-м $2ом $3ом проезде "],
            ["^ [Пп]роезд ", " проезде "],

            ["^ (\\S+н)ий [Пп]роспект ", " $1ем проспекте "],
            ["^ (\\S+)[иоы]й [Пп]роспект ", " $1ом проспекте "],
            ["^ (\\S+[иы]н) [Пп]роспект ", " $1ом проспекте "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]роспект ", " $1ом $2ем проспекте "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]роспект ", " $1ем $2ом проспекте "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]роспект ", " $1ом $2ом проспекте "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]роспект ", " $1ом $2ом проспекте "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]роспект ", " $1-м $2ем проспекте "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]роспект ", " $1-м $2ом проспекте "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]роспект ", " $1-м $2ом проспекте "],
            ["^ [Пп]роспект ", " проспекте "],

            ["^ (\\S+н)ий [Пп]утепровод ", " $1ем путепроводе "],
            ["^ (\\S+)[иоы]й [Пп]утепровод ", " $1ом путепроводе "],
            ["^ (\\S+[иы]н) [Пп]утепровод ", " $1ом путепроводе "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Пп]утепровод ", " $1ом $2ем путепроводе "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Пп]утепровод ", " $1ем $2ом путепроводе "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Пп]утепровод ", " $1ом $2ом путепроводе "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Пп]утепровод ", " $1ом $2ом путепроводе "],
            ["^ (\\d+)-й (\\S+н)ий [Пп]утепровод ", " $1-м $2ем путепроводе "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Пп]утепровод ", " $1-м $2ом путепроводе "],
            ["^ (\\d+)-й (\\S+[иы]н) [Пп]утепровод ", " $1-м $2ом путепроводе "],
            ["^ [Пп]утепровод ", " путепроводе "],

            ["^ (\\S+н)ий [Сс]пуск ", " $1ем спуске "],
            ["^ (\\S+)[иоы]й [Сс]пуск ", " $1ом спуске "],
            ["^ (\\S+[еёо]в) [Сс]пуск ", " $1ом спуске "],
            ["^ (\\S+[иы]н) [Сс]пуск ", " $1ом спуске "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]пуск ", " $1ом $2ем спуске "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]пуск ", " $1ем $2ом спуске "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]пуск ", " $1ом $2ом спуске "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Сс]пуск ", " $1ом $2ом спуске "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]пуск ", " $1ом $2ом спуске "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]пуск ", " $1-м $2ем спуске "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]пуск ", " $1-м $2ом спуске "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Сс]пуск ", " $1-м $2ом спуске "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]пуск ", " $1-м $2ом спуске "],
            ["^ [Сс]пуск ", " спуске "],

            ["^ (\\S+н)ий [Сс]ъезд ", " $1ем съезде "],
            ["^ (\\S+)[иоы]й [Сс]ъезд ", " $1ом съезде "],
            ["^ (\\S+[иы]н) [Сс]ъезд ", " $1ом съезде "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Сс]ъезд ", " $1ом $2ем съезде "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Сс]ъезд ", " $1ем $2ом съезде "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Сс]ъезд ", " $1ом $2ом съезде "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Сс]ъезд ", " $1ом $2ом съезде "],
            ["^ (\\d+)-й (\\S+н)ий [Сс]ъезд ", " $1-м $2ем съезде "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Сс]ъезд ", " $1-м $2ом съезде "],
            ["^ (\\d+)-й (\\S+[иы]н) [Сс]ъезд ", " $1-м $2ом съезде "],
            ["^ [Сс]ъезд ", " съезде "],

            ["^ (\\S+н)ий [Тт][уо]ннель ", " $1ем тоннеле "],
            ["^ (\\S+)[иоы]й [Тт][уо]ннель ", " $1ом тоннеле "],
            ["^ (\\S+[иы]н) [Тт][уо]ннель ", " $1ом тоннеле "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт][уо]ннель ", " $1ом $2ем тоннеле "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт][уо]ннель ", " $1ем $2ом тоннеле "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт][уо]ннель ", " $1ом $2ом тоннеле "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт][уо]ннель ", " $1ом $2ом тоннеле "],
            ["^ (\\d+)-й (\\S+н)ий [Тт][уо]ннель ", " $1-м $2ем тоннеле "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт][уо]ннель ", " $1-м $2ом тоннеле "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт][уо]ннель ", " $1-м $2ом тоннеле "],
            ["^ [Тт][уо]ннель ", " тоннеле "],

            ["^ (\\S+н)ий [Тт]ракт ", " $1ем тракте "],
            ["^ (\\S+)[иоы]й [Тт]ракт ", " $1ом тракте "],
            ["^ (\\S+[еёо]в) [Тт]ракт ", " $1ом тракте "],
            ["^ (\\S+[иы]н) [Тт]ракт ", " $1ом тракте "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]ракт ", " $1ом $2ем тракте "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]ракт ", " $1ем $2ом тракте "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]ракт ", " $1ом $2ом тракте "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]ракт ", " $1ом $2ом тракте "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]ракт ", " $1ом $2ом тракте "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]ракт ", " $1-м $2ем тракте "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]ракт ", " $1-м $2ом тракте "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]ракт ", " $1-м $2ом тракте "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]ракт ", " $1-м $2ом тракте "],
            ["^ [Тт]ракт ", " тракте "],

            ["^ (\\S+н)ий [Тт]упик ", " $1ем тупике "],
            ["^ (\\S+)[иоы]й [Тт]упик ", " $1ом тупике "],
            ["^ (\\S+[еёо]в) [Тт]упик ", " $1ом тупике "],
            ["^ (\\S+[иы]н) [Тт]упик ", " $1ом тупике "],
            ["^ (\\S+)[иоы]й (\\S+н)ий [Тт]упик ", " $1ом $2ем тупике "],
            ["^ (\\S+н)ий (\\S+)[иоы]й [Тт]упик ", " $1ем $2ом тупике "],
            ["^ (\\S+)[иоы]й (\\S+)[иоы]й [Тт]упик ", " $1ом $2ом тупике "],
            ["^ (\\S+)[иоы]й (\\S+[еёо]в) [Тт]упик ", " $1ом $2ом тупике "],
            ["^ (\\S+)[иоы]й (\\S+[иы]н) [Тт]упик ", " $1ом $2ом тупике "],
            ["^ (\\d+)-й [Тт]упик ", " $1-м тупике "],
            ["^ (\\d+)-й (\\S+н)ий [Тт]упик ", " $1-м $2ем тупике "],
            ["^ (\\d+)-й (\\S+)[иоы]й [Тт]упик ", " $1-м $2ом тупике "],
            ["^ (\\d+)-й (\\S+[еёо]в) [Тт]упик ", " $1-м $2ом тупике "],
            ["^ (\\d+)-й (\\S+[иы]н) [Тт]упик ", " $1-м $2ом тупике "],
            ["^ [Тт]упик ", " тупике "],

            ["^ (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1м $2кольце "],
            ["^ (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1м $2кольце "],
            ["^ (\\S+[ео])е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1м $2м $3кольце "],
            ["^ (\\S+ье) (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1м $2м $3кольце "],
            ["^ (\\d+)-е (\\S+[ео])е ([Пп]олу)?[Кк]ольцо ", " $1-м $2м $3кольце "],
            ["^ (\\d+)-е (\\S+ье) ([Пп]олу)?[Кк]ольцо ", " $1-м $2м $3кольце "],
            ["^ ([Пп]олу)?[Кк]ольцо ", " $1кольце "],

            ["^ (\\S+[ео])е [Шш]оссе ", " $1м шоссе "],
            ["^ (\\S+ье) [Шш]оссе ", " $1м шоссе "],
            ["^ (\\S+[ео])е (\\S+[ео])е [Шш]оссе ", " $1м $2м шоссе "],
            ["^ (\\S+ье) (\\S+[ео])е [Шш]оссе ", " $1м $2м шоссе "],
            ["^ (\\d+)-е (\\S+[ео])е [Шш]оссе ", " $1-м $2м шоссе "],
            ["^ (\\d+)-е (\\S+ье) [Шш]оссе ", " $1-м $2м шоссе "],

            [" ([Тт])ретом ", " $1ретьем "],
            ["([жч])ом ", "$1ьем "]
        ]
    }
}

},{}],23:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "første",
                "2": "anden",
                "3": "tredje",
                "4": "fjerde",
                "5": "femte",
                "6": "sjette",
                "7": "syvende",
                "8": "ottende",
                "9": "niende",
                "10": "tiende"
            },
            "direction": {
                "north": "Nord",
                "northeast": "Nordøst",
                "east": "Øst",
                "southeast": "Sydøst",
                "south": "Syd",
                "southwest": "Sydvest",
                "west": "Vest",
                "northwest": "Nordvest"
            },
            "modifier": {
                "left": "venstresving",
                "right": "højresving",
                "sharp left": "skarpt venstresving",
                "sharp right": "skarpt højresving",
                "slight left": "svagt venstresving",
                "slight right": "svagt højresving",
                "straight": "ligeud",
                "uturn": "U-vending"
            },
            "lanes": {
                "xo": "Hold til højre",
                "ox": "Hold til venstre",
                "xox": "Benyt midterste spor",
                "oxo": "Hold til højre eller venstre"
            }
        },
        "modes": {
            "ferry": {
                "default": "Tag færgen",
                "name": "Tag færgen {way_name}",
                "destination": "Tag færgen i retning {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} derefter, efter {distance}, {instruction_two}",
            "two linked": "{instruction_one}, derefter {instruction_two}",
            "one in distance": "Efter {distance} {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "afkørsel {exit}"
        },
        "arrive": {
            "default": {
                "default": "Du er ankommet til din {nth} destination",
                "upcoming": "Du vil ankomme til din {nth} destination",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}"
            },
            "left": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til venstre",
                "upcoming": "Du vil ankomme til din {nth} destination på venstre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til venstre"
            },
            "right": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til højre",
                "upcoming": "Du vil ankomme til din {nth} destination på højre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til højre"
            },
            "sharp left": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til venstre",
                "upcoming": "Du vil ankomme til din {nth} destination på venstre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til venstre"
            },
            "sharp right": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til højre",
                "upcoming": "Du vil ankomme til din {nth} destination på højre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til højre"
            },
            "slight right": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til højre",
                "upcoming": "Du vil ankomme til din {nth} destination på højre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til højre"
            },
            "slight left": {
                "default": "Du er ankommet til din {nth} destination, som befinder sig til venstre",
                "upcoming": "Du vil ankomme til din {nth} destination på venstre hånd",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, som befinder sig til venstre"
            },
            "straight": {
                "default": "Du er ankommet til din {nth} destination, der befinder sig lige frem",
                "upcoming": "Du vil ankomme til din {nth} destination foran dig",
                "short": "Du er ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du er ankommet til {waypoint_name}, der befinder sig lige frem"
            }
        },
        "continue": {
            "default": {
                "default": "Drej til {modifier}",
                "name": "Drej til {modifier} videre ad {way_name}",
                "destination": "Drej til {modifier} mod {destination}",
                "exit": "Drej til {modifier} ad {way_name}"
            },
            "straight": {
                "default": "Fortsæt ligeud",
                "name": "Fortsæt ligeud ad {way_name}",
                "destination": "Fortsæt mod {destination}",
                "distance": "Fortsæt {distance} ligeud",
                "namedistance": "Fortsæt {distance} ad {way_name}"
            },
            "sharp left": {
                "default": "Drej skarpt til venstre",
                "name": "Drej skarpt til venstre videre ad {way_name}",
                "destination": "Drej skarpt til venstre mod {destination}"
            },
            "sharp right": {
                "default": "Drej skarpt til højre",
                "name": "Drej skarpt til højre videre ad {way_name}",
                "destination": "Drej skarpt til højre mod {destination}"
            },
            "slight left": {
                "default": "Drej left til venstre",
                "name": "Drej let til venstre videre ad {way_name}",
                "destination": "Drej let til venstre mod {destination}"
            },
            "slight right": {
                "default": "Drej let til højre",
                "name": "Drej let til højre videre ad {way_name}",
                "destination": "Drej let til højre mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending",
                "name": "Foretag en U-vending tilbage ad {way_name}",
                "destination": "Foretag en U-vending mod {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Kør mod {direction}",
                "name": "Kør mod {direction} ad {way_name}",
                "namedistance": "Fortsæt {distance} ad {way_name}mod {direction}"
            }
        },
        "end of road": {
            "default": {
                "default": "Drej til {modifier}",
                "name": "Drej til {modifier} ad {way_name}",
                "destination": "Drej til {modifier} mof {destination}"
            },
            "straight": {
                "default": "Fortsæt ligeud",
                "name": "Fortsæt ligeud ad {way_name}",
                "destination": "Fortsæt ligeud mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending for enden af vejen",
                "name": "Foretag en U-vending ad {way_name} for enden af vejen",
                "destination": "Foretag en U-vending mod {destination} for enden af vejen"
            }
        },
        "fork": {
            "default": {
                "default": "Hold til {modifier} ved udfletningen",
                "name": "Hold mod {modifier} på {way_name}",
                "destination": "Hold mod {modifier} mod {destination}"
            },
            "slight left": {
                "default": "Hold til venstre ved udfletningen",
                "name": "Hold til venstre på {way_name}",
                "destination": "Hold til venstre mod {destination}"
            },
            "slight right": {
                "default": "Hold til højre ved udfletningen",
                "name": "Hold til højre på {way_name}",
                "destination": "Hold til højre mod {destination}"
            },
            "sharp left": {
                "default": "Drej skarpt til venstre ved udfletningen",
                "name": "Drej skarpt til venstre ad {way_name}",
                "destination": "Drej skarpt til venstre mod {destination}"
            },
            "sharp right": {
                "default": "Drej skarpt til højre ved udfletningen",
                "name": "Drej skarpt til højre ad {way_name}",
                "destination": "Drej skarpt til højre mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending",
                "name": "Foretag en U-vending ad {way_name}",
                "destination": "Foretag en U-vending mod {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Flet til {modifier}",
                "name": "Flet til {modifier} ad {way_name}",
                "destination": "Flet til {modifier} mod {destination}"
            },
            "straight": {
                "default": "Flet",
                "name": "Flet ind på {way_name}",
                "destination": "Flet ind mod {destination}"
            },
            "slight left": {
                "default": "Flet til venstre",
                "name": "Flet til venstre ad {way_name}",
                "destination": "Flet til venstre mod {destination}"
            },
            "slight right": {
                "default": "Flet til højre",
                "name": "Flet til højre ad {way_name}",
                "destination": "Flet til højre mod {destination}"
            },
            "sharp left": {
                "default": "Flet til venstre",
                "name": "Flet til venstre ad {way_name}",
                "destination": "Flet til venstre mod {destination}"
            },
            "sharp right": {
                "default": "Flet til højre",
                "name": "Flet til højre ad {way_name}",
                "destination": "Flet til højre mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending",
                "name": "Foretag en U-vending ad {way_name}",
                "destination": "Foretag en U-vending mod {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Fortsæt {modifier}",
                "name": "Fortsæt {modifier} ad {way_name}",
                "destination": "Fortsæt {modifier} mod {destination}"
            },
            "straight": {
                "default": "Fortsæt ligeud",
                "name": "Fortsæt ad {way_name}",
                "destination": "Fortsæt mod {destination}"
            },
            "sharp left": {
                "default": "Drej skarpt til venstre",
                "name": "Drej skarpt til venstre ad {way_name}",
                "destination": "Drej skarpt til venstre mod {destination}"
            },
            "sharp right": {
                "default": "Drej skarpt til højre",
                "name": "Drej skarpt til højre ad {way_name}",
                "destination": "Drej skarpt til højre mod {destination}"
            },
            "slight left": {
                "default": "Fortsæt til venstre",
                "name": "Fortsæt til venstre ad {way_name}",
                "destination": "Fortsæt til venstre mod {destination}"
            },
            "slight right": {
                "default": "Fortsæt til højre",
                "name": "Fortsæt til højre ad {way_name}",
                "destination": "Fortsæt til højre mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending",
                "name": "Foretag en U-vending ad {way_name}",
                "destination": "Foretag en U-vending mod {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Fortsæt {modifier}",
                "name": "Fortsæt {modifier} ad {way_name}",
                "destination": "Fortsæt {modifier} mod {destination}"
            },
            "uturn": {
                "default": "Foretag en U-vending",
                "name": "Foretag en U-vending ad {way_name}",
                "destination": "Foretag en U-vending mod {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Tag afkørslen",
                "name": "Tag afkørslen ad {way_name}",
                "destination": "Tag afkørslen mod {destination}",
                "exit": "Vælg afkørsel {exit}",
                "exit_destination": "Vælg afkørsel {exit} mod {destination}"
            },
            "left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}",
                "exit": "Vælg afkørsel {exit} til venstre",
                "exit_destination": "Vælg afkørsel {exit} til venstre mod {destination}\n"
            },
            "right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}",
                "exit": "Vælg afkørsel {exit} til højre",
                "exit_destination": "Vælg afkørsel {exit} til højre mod {destination}"
            },
            "sharp left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}",
                "exit": "Vælg afkørsel {exit} til venstre",
                "exit_destination": "Vælg afkørsel {exit} til venstre mod {destination}\n"
            },
            "sharp right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}",
                "exit": "Vælg afkørsel {exit} til højre",
                "exit_destination": "Vælg afkørsel {exit} til højre mod {destination}"
            },
            "slight left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}",
                "exit": "Vælg afkørsel {exit} til venstre",
                "exit_destination": "Vælg afkørsel {exit} til venstre mod {destination}\n"
            },
            "slight right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}",
                "exit": "Vælg afkørsel {exit} til højre",
                "exit_destination": "Vælg afkørsel {exit} til højre mod {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Tag afkørslen",
                "name": "Tag afkørslen ad {way_name}",
                "destination": "Tag afkørslen mod {destination}"
            },
            "left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}"
            },
            "right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}"
            },
            "sharp left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}"
            },
            "sharp right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}"
            },
            "slight left": {
                "default": "Tag afkørslen til venstre",
                "name": "Tag afkørslen til venstre ad {way_name}",
                "destination": "Tag afkørslen til venstre mod {destination}"
            },
            "slight right": {
                "default": "Tag afkørslen til højre",
                "name": "Tag afkørslen til højre ad {way_name}",
                "destination": "Tag afkørslen til højre mod {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Kør ind i rundkørslen",
                    "name": "Tag rundkørslen og kør fra ad {way_name}",
                    "destination": "Tag rundkørslen og kør mod {destination}"
                },
                "name": {
                    "default": "Kør ind i {rotary_name}",
                    "name": "Kør ind i {rotary_name} og kør ad {way_name} ",
                    "destination": "Kør ind i {rotary_name} og kør mod {destination}"
                },
                "exit": {
                    "default": "Tag rundkørslen og forlad ved {exit_number} afkørsel",
                    "name": "Tag rundkørslen og forlad ved {exit_number} afkørsel ad {way_name}",
                    "destination": "Tag rundkørslen og forlad ved {exit_number} afkørsel mod {destination}"
                },
                "name_exit": {
                    "default": "Kør ind i {rotary_name} og forlad ved {exit_number} afkørsel",
                    "name": "Kør ind i {rotary_name} og forlad ved {exit_number} afkørsel ad {way_name}",
                    "destination": "Kør ind i {rotary_name} og forlad ved {exit_number} afkørsel mod {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Tag rundkørslen og forlad ved {exit_number} afkørsel",
                    "name": "Tag rundkørslen og forlad ved {exit_number} afkørsel ad {way_name}",
                    "destination": "Tag rundkørslen og forlad ved {exit_number} afkørsel mod {destination}"
                },
                "default": {
                    "default": "Kør ind i rundkørslen",
                    "name": "Tag rundkørslen og kør fra ad {way_name}",
                    "destination": "Tag rundkørslen og kør mod {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Foretag et {modifier}",
                "name": "Foretag et {modifier} ad {way_name}",
                "destination": "Foretag et {modifier} mod {destination}"
            },
            "left": {
                "default": "Drej til venstre",
                "name": "Drej til venstre ad {way_name}",
                "destination": "Drej til venstre mod {destination}"
            },
            "right": {
                "default": "Drej til højre",
                "name": "Drej til højre ad {way_name}",
                "destination": "Drej til højre mod {destination}"
            },
            "straight": {
                "default": "Fortsæt ligeud",
                "name": "Fortsæt ligeud ad {way_name}",
                "destination": "Fortsæt ligeud mod {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Forlad rundkørslen",
                "name": "Forlad rundkørslen ad {way_name}",
                "destination": "Forlad rundkørslen mod  {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Forlad rundkørslen",
                "name": "Forlad rundkørslen ad {way_name}",
                "destination": "Forlad rundkørslen mod {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Foretag et {modifier}",
                "name": "Foretag et {modifier} ad {way_name}",
                "destination": "Foretag et {modifier} mod {destination}"
            },
            "left": {
                "default": "Drej til venstre",
                "name": "Drej til venstre ad {way_name}",
                "destination": "Drej til venstre mod {destination}"
            },
            "right": {
                "default": "Drej til højre",
                "name": "Drej til højre ad {way_name}",
                "destination": "Drej til højre mod {destination}"
            },
            "straight": {
                "default": "Fortsæt ligeud",
                "name": "Kør ligeud ad {way_name}",
                "destination": "Kør ligeud mod {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Fortsæt ligeud"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],24:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "erste",
                "2": "zweite",
                "3": "dritte",
                "4": "vierte",
                "5": "fünfte",
                "6": "sechste",
                "7": "siebente",
                "8": "achte",
                "9": "neunte",
                "10": "zehnte"
            },
            "direction": {
                "north": "Norden",
                "northeast": "Nordosten",
                "east": "Osten",
                "southeast": "Südosten",
                "south": "Süden",
                "southwest": "Südwesten",
                "west": "Westen",
                "northwest": "Nordwesten"
            },
            "modifier": {
                "left": "links",
                "right": "rechts",
                "sharp left": "scharf links",
                "sharp right": "scharf rechts",
                "slight left": "leicht links",
                "slight right": "leicht rechts",
                "straight": "geradeaus",
                "uturn": "180°-Wendung"
            },
            "lanes": {
                "xo": "Rechts halten",
                "ox": "Links halten",
                "xox": "Mittlere Spur nutzen",
                "oxo": "Rechts oder links halten"
            }
        },
        "modes": {
            "ferry": {
                "default": "Fähre nehmen",
                "name": "Fähre nehmen {way_name}",
                "destination": "Fähre nehmen Richtung {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} danach in {distance} {instruction_two}",
            "two linked": "{instruction_one} danach {instruction_two}",
            "one in distance": "In {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "Sie haben Ihr {nth} Ziel erreicht",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}"
            },
            "left": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich links"
            },
            "right": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich rechts"
            },
            "sharp left": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich links"
            },
            "sharp right": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich rechts"
            },
            "slight right": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich rechts",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich rechts"
            },
            "slight left": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich links",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich links"
            },
            "straight": {
                "default": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich geradeaus",
                "upcoming": "Sie haben Ihr {nth} Ziel erreicht, es befindet sich geradeaus",
                "short": "Sie haben Ihr {nth} Ziel erreicht",
                "short-upcoming": "Sie haben Ihr {nth} Ziel erreicht",
                "named": "Sie haben Ihr {waypoint_name}, es befindet sich geradeaus"
            }
        },
        "continue": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} weiterfahren auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}",
                "exit": "{modifier} abbiegen auf {way_name}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Weiterfahren in Richtung {destination}",
                "distance": "Geradeaus weiterfahren für {distance}",
                "namedistance": "Geradeaus weiterfahren auf {way_name} für {distance}"
            },
            "sharp left": {
                "default": "Scharf links",
                "name": "Scharf links weiterfahren auf {way_name}",
                "destination": "Scharf links Richtung {destination}"
            },
            "sharp right": {
                "default": "Scharf rechts",
                "name": "Scharf rechts weiterfahren auf {way_name}",
                "destination": "Scharf rechts Richtung {destination}"
            },
            "slight left": {
                "default": "Leicht links",
                "name": "Leicht links weiter auf {way_name}",
                "destination": "Leicht links weiter Richtung {destination}"
            },
            "slight right": {
                "default": "Leicht rechts weiter",
                "name": "Leicht rechts weiter auf {way_name}",
                "destination": "Leicht rechts weiter Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung",
                "name": "180°-Wendung auf {way_name}",
                "destination": "180°-Wendung Richtung {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Fahren Sie Richtung {direction}",
                "name": "Fahren Sie Richtung {direction} auf {way_name}",
                "namedistance": "Fahren Sie Richtung {direction} auf {way_name} für {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} abbiegen auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Geradeaus weiterfahren Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung am Ende der Straße",
                "name": "180°-Wendung auf {way_name} am Ende der Straße",
                "destination": "180°-Wendung Richtung {destination} am Ende der Straße"
            }
        },
        "fork": {
            "default": {
                "default": "{modifier} halten an der Gabelung",
                "name": "{modifier} halten an der Gabelung auf {way_name}",
                "destination": "{modifier}  halten an der Gabelung Richtung {destination}"
            },
            "slight left": {
                "default": "Links halten an der Gabelung",
                "name": "Links halten an der Gabelung auf {way_name}",
                "destination": "Links halten an der Gabelung Richtung {destination}"
            },
            "slight right": {
                "default": "Rechts halten an der Gabelung",
                "name": "Rechts halten an der Gabelung auf {way_name}",
                "destination": "Rechts halten an der Gabelung Richtung {destination}"
            },
            "sharp left": {
                "default": "Scharf links abbiegen an der Gabelung",
                "name": "Scharf links auf {way_name}",
                "destination": "Scharf links Richtung {destination}"
            },
            "sharp right": {
                "default": "Scharf rechts abbiegen an der Gabelung",
                "name": "Scharf rechts auf {way_name}",
                "destination": "Scharf rechts Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung",
                "name": "180°-Wendung auf {way_name}",
                "destination": "180°-Wendung Richtung {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "{modifier} auffahren",
                "name": "{modifier} auffahren auf {way_name}",
                "destination": "{modifier} auffahren Richtung {destination}"
            },
            "straight": {
                "default": "geradeaus auffahren",
                "name": "geradeaus auffahren auf {way_name}",
                "destination": "geradeaus auffahren Richtung {destination}"
            },
            "slight left": {
                "default": "Leicht links auffahren",
                "name": "Leicht links auffahren auf {way_name}",
                "destination": "Leicht links auffahren Richtung {destination}"
            },
            "slight right": {
                "default": "Leicht rechts auffahren",
                "name": "Leicht rechts auffahren auf {way_name}",
                "destination": "Leicht rechts auffahren Richtung {destination}"
            },
            "sharp left": {
                "default": "Scharf links auffahren",
                "name": "Scharf links auffahren auf {way_name}",
                "destination": "Scharf links auffahren Richtung {destination}"
            },
            "sharp right": {
                "default": "Scharf rechts auffahren",
                "name": "Scharf rechts auffahren auf {way_name}",
                "destination": "Scharf rechts auffahren Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung",
                "name": "180°-Wendung auf {way_name}",
                "destination": "180°-Wendung Richtung {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "{modifier} weiterfahren",
                "name": "{modifier} weiterfahren auf {way_name}",
                "destination": "{modifier} weiterfahren Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Weiterfahren auf {way_name}",
                "destination": "Weiterfahren in Richtung {destination}"
            },
            "sharp left": {
                "default": "Scharf links",
                "name": "Scharf links auf {way_name}",
                "destination": "Scharf links Richtung {destination}"
            },
            "sharp right": {
                "default": "Scharf rechts",
                "name": "Scharf rechts auf {way_name}",
                "destination": "Scharf rechts Richtung {destination}"
            },
            "slight left": {
                "default": "Leicht links weiter",
                "name": "Leicht links weiter auf {way_name}",
                "destination": "Leicht links weiter Richtung {destination}"
            },
            "slight right": {
                "default": "Leicht rechts weiter",
                "name": "Leicht rechts weiter auf {way_name}",
                "destination": "Leicht rechts weiter Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung",
                "name": "180°-Wendung auf {way_name}",
                "destination": "180°-Wendung Richtung {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "{modifier} weiterfahren",
                "name": "{modifier} weiterfahren auf {way_name}",
                "destination": "{modifier} weiterfahren Richtung {destination}"
            },
            "uturn": {
                "default": "180°-Wendung",
                "name": "180°-Wendung auf {way_name}",
                "destination": "180°-Wendung Richtung {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Ausfahrt nehmen",
                "name": "Ausfahrt nehmen auf {way_name}",
                "destination": "Ausfahrt nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} nehmen",
                "exit_destination": "Ausfahrt {exit} nehmen Richtung {destination}"
            },
            "left": {
                "default": "Ausfahrt links nehmen",
                "name": "Ausfahrt links nehmen auf {way_name}",
                "destination": "Ausfahrt links nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} links nehmen",
                "exit_destination": "Ausfahrt {exit} links nehmen Richtung {destination}"
            },
            "right": {
                "default": "Ausfahrt rechts nehmen",
                "name": "Ausfahrt rechts nehmen Richtung {way_name}",
                "destination": "Ausfahrt rechts nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} rechts nehmen",
                "exit_destination": "Ausfahrt {exit} nehmen Richtung {destination}"
            },
            "sharp left": {
                "default": "Ausfahrt links nehmen",
                "name": "Ausfahrt links Seite nehmen auf {way_name}",
                "destination": "Ausfahrt links nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} links nehmen",
                "exit_destination": "Ausfahrt{exit} links nehmen Richtung {destination}"
            },
            "sharp right": {
                "default": "Ausfahrt rechts nehmen",
                "name": "Ausfahrt rechts nehmen auf {way_name}",
                "destination": "Ausfahrt rechts nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} rechts nehmen",
                "exit_destination": "Ausfahrt {exit} nehmen Richtung {destination}"
            },
            "slight left": {
                "default": "Ausfahrt links nehmen",
                "name": "Ausfahrt links nehmen auf {way_name}",
                "destination": "Ausfahrt links nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} nehmen",
                "exit_destination": "Ausfahrt {exit} links nehmen Richtung {destination}"
            },
            "slight right": {
                "default": "Ausfahrt rechts nehmen",
                "name": "Ausfahrt rechts nehmen auf {way_name}",
                "destination": "Ausfahrt rechts nehmen Richtung {destination}",
                "exit": "Ausfahrt {exit} rechts nehmen",
                "exit_destination": "Ausfahrt {exit} nehmen Richtung {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Auffahrt nehmen",
                "name": "Auffahrt nehmen auf {way_name}",
                "destination": "Auffahrt nehmen Richtung {destination}"
            },
            "left": {
                "default": "Auffahrt links nehmen",
                "name": "Auffahrt links nehmen auf {way_name}",
                "destination": "Auffahrt links nehmen Richtung {destination}"
            },
            "right": {
                "default": "Auffahrt rechts nehmen",
                "name": "Auffahrt rechts nehmen auf {way_name}",
                "destination": "Auffahrt rechts nehmen Richtung {destination}"
            },
            "sharp left": {
                "default": "Auffahrt links nehmen",
                "name": "Auffahrt links nehmen auf {way_name}",
                "destination": "Auffahrt links nehmen Richtung {destination}"
            },
            "sharp right": {
                "default": "Auffahrt rechts nehmen",
                "name": "Auffahrt rechts nehmen auf {way_name}",
                "destination": "Auffahrt rechts nehmen Richtung {destination}"
            },
            "slight left": {
                "default": "Auffahrt links Seite nehmen",
                "name": "Auffahrt links nehmen auf {way_name}",
                "destination": "Auffahrt links nehmen Richtung {destination}"
            },
            "slight right": {
                "default": "Auffahrt rechts nehmen",
                "name": "Auffahrt rechts nehmen auf {way_name}",
                "destination": "Auffahrt rechts nehmen Richtung {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "In den Kreisverkehr fahren",
                    "name": "Im Kreisverkehr die Ausfahrt auf {way_name} nehmen",
                    "destination": "Im Kreisverkehr die Ausfahrt Richtung {destination} nehmen"
                },
                "name": {
                    "default": "In {rotary_name} fahren",
                    "name": "In {rotary_name} die Ausfahrt auf {way_name} nehmen",
                    "destination": "In {rotary_name} die Ausfahrt Richtung {destination} nehmen"
                },
                "exit": {
                    "default": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen",
                    "name": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen auf {way_name}",
                    "destination": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen Richtung {destination}"
                },
                "name_exit": {
                    "default": "In den Kreisverkehr fahren und {exit_number} Ausfahrt nehmen",
                    "name": "In den Kreisverkehr fahren und {exit_number} Ausfahrt nehmen auf {way_name}",
                    "destination": "In den Kreisverkehr fahren und {exit_number} Ausfahrt nehmen Richtung {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen",
                    "name": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen auf {way_name}",
                    "destination": "Im Kreisverkehr die {exit_number} Ausfahrt nehmen Richtung {destination}"
                },
                "default": {
                    "default": "In den Kreisverkehr fahren",
                    "name": "Im Kreisverkehr die Ausfahrt auf {way_name} nehmen",
                    "destination": "Im Kreisverkehr die Ausfahrt Richtung {destination} nehmen"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} abbiegen auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}"
            },
            "left": {
                "default": "Links abbiegen",
                "name": "Links abbiegen auf {way_name}",
                "destination": "Links abbiegen Richtung {destination}"
            },
            "right": {
                "default": "Rechts abbiegen",
                "name": "Rechts abbiegen auf {way_name}",
                "destination": "Rechts abbiegen Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Geradeaus weiterfahren Richtung {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} abbiegen auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}"
            },
            "left": {
                "default": "Links abbiegen",
                "name": "Links abbiegen auf {way_name}",
                "destination": "Links abbiegen Richtung {destination}"
            },
            "right": {
                "default": "Rechts abbiegen",
                "name": "Rechts abbiegen auf {way_name}",
                "destination": "Rechts abbiegen Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Geradeaus weiterfahren Richtung {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} abbiegen auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}"
            },
            "left": {
                "default": "Links abbiegen",
                "name": "Links abbiegen auf {way_name}",
                "destination": "Links abbiegen Richtung {destination}"
            },
            "right": {
                "default": "Rechts abbiegen",
                "name": "Rechts abbiegen auf {way_name}",
                "destination": "Rechts abbiegen Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Geradeaus weiterfahren Richtung {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "{modifier} abbiegen",
                "name": "{modifier} abbiegen auf {way_name}",
                "destination": "{modifier} abbiegen Richtung {destination}"
            },
            "left": {
                "default": "Links abbiegen",
                "name": "Links abbiegen auf {way_name}",
                "destination": "Links abbiegen Richtung {destination}"
            },
            "right": {
                "default": "Rechts abbiegen",
                "name": "Rechts abbiegen auf {way_name}",
                "destination": "Rechts abbiegen Richtung {destination}"
            },
            "straight": {
                "default": "Geradeaus weiterfahren",
                "name": "Geradeaus weiterfahren auf {way_name}",
                "destination": "Geradeaus weiterfahren Richtung {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Geradeaus weiterfahren"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],25:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1st",
                "2": "2nd",
                "3": "3rd",
                "4": "4th",
                "5": "5th",
                "6": "6th",
                "7": "7th",
                "8": "8th",
                "9": "9th",
                "10": "10th"
            },
            "direction": {
                "north": "north",
                "northeast": "northeast",
                "east": "east",
                "southeast": "southeast",
                "south": "south",
                "southwest": "southwest",
                "west": "west",
                "northwest": "northwest"
            },
            "modifier": {
                "left": "left",
                "right": "right",
                "sharp left": "sharp left",
                "sharp right": "sharp right",
                "slight left": "slight left",
                "slight right": "slight right",
                "straight": "straight",
                "uturn": "U-turn"
            },
            "lanes": {
                "xo": "Keep right",
                "ox": "Keep left",
                "xox": "Keep in the middle",
                "oxo": "Keep left or right"
            }
        },
        "modes": {
            "ferry": {
                "default": "Take the ferry",
                "name": "Take the ferry {way_name}",
                "destination": "Take the ferry towards {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, then, in {distance}, {instruction_two}",
            "two linked": "{instruction_one}, then {instruction_two}",
            "one in distance": "In {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "You have arrived at your {nth} destination",
                "upcoming": "You will arrive at your {nth} destination",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}"
            },
            "left": {
                "default": "You have arrived at your {nth} destination, on the left",
                "upcoming": "You will arrive at your {nth} destination, on the left",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the left"
            },
            "right": {
                "default": "You have arrived at your {nth} destination, on the right",
                "upcoming": "You will arrive at your {nth} destination, on the right",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the right"
            },
            "sharp left": {
                "default": "You have arrived at your {nth} destination, on the left",
                "upcoming": "You will arrive at your {nth} destination, on the left",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the left"
            },
            "sharp right": {
                "default": "You have arrived at your {nth} destination, on the right",
                "upcoming": "You will arrive at your {nth} destination, on the right",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the right"
            },
            "slight right": {
                "default": "You have arrived at your {nth} destination, on the right",
                "upcoming": "You will arrive at your {nth} destination, on the right",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the right"
            },
            "slight left": {
                "default": "You have arrived at your {nth} destination, on the left",
                "upcoming": "You will arrive at your {nth} destination, on the left",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, on the left"
            },
            "straight": {
                "default": "You have arrived at your {nth} destination, straight ahead",
                "upcoming": "You will arrive at your {nth} destination, straight ahead",
                "short": "You have arrived",
                "short-upcoming": "You will arrive",
                "named": "You have arrived at {waypoint_name}, straight ahead"
            }
        },
        "continue": {
            "default": {
                "default": "Turn {modifier}",
                "name": "Turn {modifier} to stay on {way_name}",
                "destination": "Turn {modifier} towards {destination}",
                "exit": "Turn {modifier} onto {way_name}"
            },
            "straight": {
                "default": "Continue straight",
                "name": "Continue straight to stay on {way_name}",
                "destination": "Continue towards {destination}",
                "distance": "Continue straight for {distance}",
                "namedistance": "Continue on {way_name} for {distance}"
            },
            "sharp left": {
                "default": "Make a sharp left",
                "name": "Make a sharp left to stay on {way_name}",
                "destination": "Make a sharp left towards {destination}"
            },
            "sharp right": {
                "default": "Make a sharp right",
                "name": "Make a sharp right to stay on {way_name}",
                "destination": "Make a sharp right towards {destination}"
            },
            "slight left": {
                "default": "Make a slight left",
                "name": "Make a slight left to stay on {way_name}",
                "destination": "Make a slight left towards {destination}"
            },
            "slight right": {
                "default": "Make a slight right",
                "name": "Make a slight right to stay on {way_name}",
                "destination": "Make a slight right towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn",
                "name": "Make a U-turn and continue on {way_name}",
                "destination": "Make a U-turn towards {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Head {direction}",
                "name": "Head {direction} on {way_name}",
                "namedistance": "Head {direction} on {way_name} for {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Turn {modifier}",
                "name": "Turn {modifier} onto {way_name}",
                "destination": "Turn {modifier} towards {destination}"
            },
            "straight": {
                "default": "Continue straight",
                "name": "Continue straight onto {way_name}",
                "destination": "Continue straight towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn at the end of the road",
                "name": "Make a U-turn onto {way_name} at the end of the road",
                "destination": "Make a U-turn towards {destination} at the end of the road"
            }
        },
        "fork": {
            "default": {
                "default": "Keep {modifier} at the fork",
                "name": "Keep {modifier} onto {way_name}",
                "destination": "Keep {modifier} towards {destination}"
            },
            "slight left": {
                "default": "Keep left at the fork",
                "name": "Keep left onto {way_name}",
                "destination": "Keep left towards {destination}"
            },
            "slight right": {
                "default": "Keep right at the fork",
                "name": "Keep right onto {way_name}",
                "destination": "Keep right towards {destination}"
            },
            "sharp left": {
                "default": "Take a sharp left at the fork",
                "name": "Take a sharp left onto {way_name}",
                "destination": "Take a sharp left towards {destination}"
            },
            "sharp right": {
                "default": "Take a sharp right at the fork",
                "name": "Take a sharp right onto {way_name}",
                "destination": "Take a sharp right towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn",
                "name": "Make a U-turn onto {way_name}",
                "destination": "Make a U-turn towards {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Merge {modifier}",
                "name": "Merge {modifier} onto {way_name}",
                "destination": "Merge {modifier} towards {destination}"
            },
            "straight": {
                "default": "Merge",
                "name": "Merge onto {way_name}",
                "destination": "Merge towards {destination}"
            },
            "slight left": {
                "default": "Merge left",
                "name": "Merge left onto {way_name}",
                "destination": "Merge left towards {destination}"
            },
            "slight right": {
                "default": "Merge right",
                "name": "Merge right onto {way_name}",
                "destination": "Merge right towards {destination}"
            },
            "sharp left": {
                "default": "Merge left",
                "name": "Merge left onto {way_name}",
                "destination": "Merge left towards {destination}"
            },
            "sharp right": {
                "default": "Merge right",
                "name": "Merge right onto {way_name}",
                "destination": "Merge right towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn",
                "name": "Make a U-turn onto {way_name}",
                "destination": "Make a U-turn towards {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} onto {way_name}",
                "destination": "Continue {modifier} towards {destination}"
            },
            "straight": {
                "default": "Continue straight",
                "name": "Continue onto {way_name}",
                "destination": "Continue towards {destination}"
            },
            "sharp left": {
                "default": "Take a sharp left",
                "name": "Take a sharp left onto {way_name}",
                "destination": "Take a sharp left towards {destination}"
            },
            "sharp right": {
                "default": "Take a sharp right",
                "name": "Take a sharp right onto {way_name}",
                "destination": "Take a sharp right towards {destination}"
            },
            "slight left": {
                "default": "Continue slightly left",
                "name": "Continue slightly left onto {way_name}",
                "destination": "Continue slightly left towards {destination}"
            },
            "slight right": {
                "default": "Continue slightly right",
                "name": "Continue slightly right onto {way_name}",
                "destination": "Continue slightly right towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn",
                "name": "Make a U-turn onto {way_name}",
                "destination": "Make a U-turn towards {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} onto {way_name}",
                "destination": "Continue {modifier} towards {destination}"
            },
            "uturn": {
                "default": "Make a U-turn",
                "name": "Make a U-turn onto {way_name}",
                "destination": "Make a U-turn towards {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Take the ramp",
                "name": "Take the ramp onto {way_name}",
                "destination": "Take the ramp towards {destination}",
                "exit": "Take exit {exit}",
                "exit_destination": "Take exit {exit} towards {destination}"
            },
            "left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            },
            "sharp left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "sharp right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            },
            "slight left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "slight right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Take the ramp",
                "name": "Take the ramp onto {way_name}",
                "destination": "Take the ramp towards {destination}"
            },
            "left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}"
            },
            "right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}"
            },
            "sharp left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}"
            },
            "sharp right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}"
            },
            "slight left": {
                "default": "Take the ramp on the left",
                "name": "Take the ramp on the left onto {way_name}",
                "destination": "Take the ramp on the left towards {destination}"
            },
            "slight right": {
                "default": "Take the ramp on the right",
                "name": "Take the ramp on the right onto {way_name}",
                "destination": "Take the ramp on the right towards {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Enter the traffic circle",
                    "name": "Enter the traffic circle and exit onto {way_name}",
                    "destination": "Enter the traffic circle and exit towards {destination}"
                },
                "name": {
                    "default": "Enter {rotary_name}",
                    "name": "Enter {rotary_name} and exit onto {way_name}",
                    "destination": "Enter {rotary_name} and exit towards {destination}"
                },
                "exit": {
                    "default": "Enter the traffic circle and take the {exit_number} exit",
                    "name": "Enter the traffic circle and take the {exit_number} exit onto {way_name}",
                    "destination": "Enter the traffic circle and take the {exit_number} exit towards {destination}"
                },
                "name_exit": {
                    "default": "Enter {rotary_name} and take the {exit_number} exit",
                    "name": "Enter {rotary_name} and take the {exit_number} exit onto {way_name}",
                    "destination": "Enter {rotary_name} and take the {exit_number} exit towards {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Enter the traffic circle and take the {exit_number} exit",
                    "name": "Enter the traffic circle and take the {exit_number} exit onto {way_name}",
                    "destination": "Enter the traffic circle and take the {exit_number} exit towards {destination}"
                },
                "default": {
                    "default": "Enter the traffic circle",
                    "name": "Enter the traffic circle and exit onto {way_name}",
                    "destination": "Enter the traffic circle and exit towards {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Make a {modifier}",
                "name": "Make a {modifier} onto {way_name}",
                "destination": "Make a {modifier} towards {destination}"
            },
            "left": {
                "default": "Turn left",
                "name": "Turn left onto {way_name}",
                "destination": "Turn left towards {destination}"
            },
            "right": {
                "default": "Turn right",
                "name": "Turn right onto {way_name}",
                "destination": "Turn right towards {destination}"
            },
            "straight": {
                "default": "Continue straight",
                "name": "Continue straight onto {way_name}",
                "destination": "Continue straight towards {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Exit the traffic circle",
                "name": "Exit the traffic circle onto {way_name}",
                "destination": "Exit the traffic circle towards {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Exit the traffic circle",
                "name": "Exit the traffic circle onto {way_name}",
                "destination": "Exit the traffic circle towards {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Make a {modifier}",
                "name": "Make a {modifier} onto {way_name}",
                "destination": "Make a {modifier} towards {destination}"
            },
            "left": {
                "default": "Turn left",
                "name": "Turn left onto {way_name}",
                "destination": "Turn left towards {destination}"
            },
            "right": {
                "default": "Turn right",
                "name": "Turn right onto {way_name}",
                "destination": "Turn right towards {destination}"
            },
            "straight": {
                "default": "Go straight",
                "name": "Go straight onto {way_name}",
                "destination": "Go straight towards {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continue straight"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],26:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1.",
                "2": "2.",
                "3": "3.",
                "4": "4.",
                "5": "5.",
                "6": "6.",
                "7": "7.",
                "8": "8.",
                "9": "9.",
                "10": "10."
            },
            "direction": {
                "north": "norden",
                "northeast": "nord-orienten",
                "east": "orienten",
                "southeast": "sud-orienten",
                "south": "suden",
                "southwest": "sud-okcidenten",
                "west": "okcidenten",
                "northwest": "nord-okcidenten"
            },
            "modifier": {
                "left": "maldekstren",
                "right": "dekstren",
                "sharp left": "maldekstregen",
                "sharp right": "dekstregen",
                "slight left": "maldekstreten",
                "slight right": "dekstreten",
                "straight": "rekten",
                "uturn": "turniĝu malantaŭen"
            },
            "lanes": {
                "xo": "Veturu dekstre",
                "ox": "Veturu maldekstre",
                "xox": "Veturu meze",
                "oxo": "Veturu dekstre aŭ maldekstre"
            }
        },
        "modes": {
            "ferry": {
                "default": "Enpramiĝu",
                "name": "Enpramiĝu {way_name}",
                "destination": "Enpramiĝu direkte al {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} kaj post {distance} {instruction_two}",
            "two linked": "{instruction_one} kaj sekve {instruction_two}",
            "one in distance": "Post {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "elveturejo {exit}"
        },
        "arrive": {
            "default": {
                "default": "Vi atingis vian {nth} celon",
                "upcoming": "Vi atingos vian {nth} celon",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}"
            },
            "left": {
                "default": "Vi atingis vian {nth} celon ĉe maldekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe maldekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe maldekstre"
            },
            "right": {
                "default": "Vi atingis vian {nth} celon ĉe dekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe dekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe dekstre"
            },
            "sharp left": {
                "default": "Vi atingis vian {nth} celon ĉe maldekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe maldekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe maldekstre"
            },
            "sharp right": {
                "default": "Vi atingis vian {nth} celon ĉe dekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe dekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe dekstre"
            },
            "slight right": {
                "default": "Vi atingis vian {nth} celon ĉe dekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe dekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe dekstre"
            },
            "slight left": {
                "default": "Vi atingis vian {nth} celon ĉe maldekstre",
                "upcoming": "Vi atingos vian {nth} celon ĉe maldekstre",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name}, ĉe maldekstre"
            },
            "straight": {
                "default": "Vi atingis vian {nth} celon",
                "upcoming": "Vi atingos vian {nth} celon rekte",
                "short": "Vi atingis",
                "short-upcoming": "Vi atingos",
                "named": "Vi atingis {waypoint_name} antaŭe"
            }
        },
        "continue": {
            "default": {
                "default": "Veturu {modifier}",
                "name": "Veturu {modifier} al {way_name}",
                "destination": "Veturu {modifier} direkte al {destination}",
                "exit": "Veturu {modifier} direkte al {way_name}"
            },
            "straight": {
                "default": "Veturu rekten",
                "name": "Veturu rekten al {way_name}",
                "destination": "Veturu rekten direkte al {destination}",
                "distance": "Veturu rekten dum {distance}",
                "namedistance": "Veturu rekten al {way_name} dum {distance}"
            },
            "sharp left": {
                "default": "Turniĝu ege maldekstren",
                "name": "Turniĝu ege maldekstren al {way_name}",
                "destination": "Turniĝu ege maldekstren direkte al {destination}"
            },
            "sharp right": {
                "default": "Turniĝu ege dekstren",
                "name": "Turniĝu ege dekstren al {way_name}",
                "destination": "Turniĝu ege dekstren direkte al {destination}"
            },
            "slight left": {
                "default": "Turniĝu ete maldekstren",
                "name": "Turniĝu ete maldekstren al {way_name}",
                "destination": "Turniĝu ete maldekstren direkte al {destination}"
            },
            "slight right": {
                "default": "Turniĝu ete dekstren",
                "name": "Turniĝu ete dekstren al {way_name}",
                "destination": "Turniĝu ete dekstren direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen",
                "name": "Turniĝu malantaŭen al {way_name}",
                "destination": "Turniĝu malantaŭen direkte al {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Direktiĝu {direction}",
                "name": "Direktiĝu {direction} al {way_name}",
                "namedistance": "Direktiĝu {direction} al {way_name} tra {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Veturu {modifier}",
                "name": "Veturu {modifier} direkte al {way_name}",
                "destination": "Veturu {modifier} direkte al {destination}"
            },
            "straight": {
                "default": "Veturu rekten",
                "name": "Veturu rekten al {way_name}",
                "destination": "Veturu rekten direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen ĉe fino de la vojo",
                "name": "Turniĝu malantaŭen al {way_name} ĉe fino de la vojo",
                "destination": "Turniĝu malantaŭen direkte al {destination} ĉe fino de la vojo"
            }
        },
        "fork": {
            "default": {
                "default": "Daŭru {modifier} ĉe la vojforko",
                "name": "Pluu {modifier} al {way_name}",
                "destination": "Pluu {modifier} direkte al {destination}"
            },
            "slight left": {
                "default": "Maldekstren ĉe la vojforko",
                "name": "Pluu maldekstren al {way_name}",
                "destination": "Pluu maldekstren direkte al {destination}"
            },
            "slight right": {
                "default": "Dekstren ĉe la vojforko",
                "name": "Pluu dekstren al {way_name}",
                "destination": "Pluu dekstren direkte al {destination}"
            },
            "sharp left": {
                "default": "Ege maldekstren ĉe la vojforko",
                "name": "Turniĝu ege maldekstren al {way_name}",
                "destination": "Turniĝu ege maldekstren direkte al {destination}"
            },
            "sharp right": {
                "default": "Ege dekstren ĉe la vojforko",
                "name": "Turniĝu ege dekstren al {way_name}",
                "destination": "Turniĝu ege dekstren direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen",
                "name": "Turniĝu malantaŭen al {way_name}",
                "destination": "Turniĝu malantaŭen direkte al {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Enveturu {modifier}",
                "name": "Enveturu {modifier} al {way_name}",
                "destination": "Enveturu {modifier} direkte al {destination}"
            },
            "straight": {
                "default": "Enveturu",
                "name": "Enveturu al {way_name}",
                "destination": "Enveturu direkte al {destination}"
            },
            "slight left": {
                "default": "Enveturu de maldekstre",
                "name": "Enveturu de maldekstre al {way_name}",
                "destination": "Enveturu de maldekstre direkte al {destination}"
            },
            "slight right": {
                "default": "Enveturu de dekstre",
                "name": "Enveturu de dekstre al {way_name}",
                "destination": "Enveturu de dekstre direkte al {destination}"
            },
            "sharp left": {
                "default": "Enveturu de maldekstre",
                "name": "Enveture de maldekstre al {way_name}",
                "destination": "Enveturu de maldekstre direkte al {destination}"
            },
            "sharp right": {
                "default": "Enveturu de dekstre",
                "name": "Enveturu de dekstre al {way_name}",
                "destination": "Enveturu de dekstre direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen",
                "name": "Turniĝu malantaŭen al {way_name}",
                "destination": "Turniĝu malantaŭen direkte al {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Pluu {modifier}",
                "name": "Pluu {modifier} al {way_name}",
                "destination": "Pluu {modifier} direkte al {destination}"
            },
            "straight": {
                "default": "Veturu rekten",
                "name": "Veturu rekten al {way_name}",
                "destination": "Veturu rekten direkte al {destination}"
            },
            "sharp left": {
                "default": "Turniĝu ege maldekstren",
                "name": "Turniĝu ege maldekstren al {way_name}",
                "destination": "Turniĝu ege maldekstren direkte al {destination}"
            },
            "sharp right": {
                "default": "Turniĝu ege dekstren",
                "name": "Turniĝu ege dekstren al {way_name}",
                "destination": "Turniĝu ege dekstren direkte al {destination}"
            },
            "slight left": {
                "default": "Pluu ete maldekstren",
                "name": "Pluu ete maldekstren al {way_name}",
                "destination": "Pluu ete maldekstren direkte al {destination}"
            },
            "slight right": {
                "default": "Pluu ete dekstren",
                "name": "Pluu ete dekstren al {way_name}",
                "destination": "Pluu ete dekstren direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen",
                "name": "Turniĝu malantaŭen al {way_name}",
                "destination": "Turniĝu malantaŭen direkte al {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Pluu {modifier}",
                "name": "Pluu {modifier} al {way_name}",
                "destination": "Pluu {modifier} direkte al {destination}"
            },
            "uturn": {
                "default": "Turniĝu malantaŭen",
                "name": "Turniĝu malantaŭen al {way_name}",
                "destination": "Turniĝu malantaŭen direkte al {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Direktiĝu al enveturejo",
                "name": "Direktiĝu al enveturejo al {way_name}",
                "destination": "Direktiĝu al enveturejo direkte al {destination}",
                "exit": "Direktiĝu al elveturejo {exit}",
                "exit_destination": "Direktiĝu al elveturejo {exit} direkte al {destination}"
            },
            "left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}",
                "exit": "Direktiĝu al elveturejo {exit} ĉe maldekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe maldekstre direkte al {destination}"
            },
            "right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}",
                "exit": "Direktiĝu al {exit} elveturejo ĉe ldekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe dekstre direkte al {destination}"
            },
            "sharp left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}",
                "exit": "Direktiĝu al {exit} elveturejo ĉe maldekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe maldekstre direkte al {destination}"
            },
            "sharp right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}",
                "exit": "Direktiĝu al elveturejo {exit} ĉe dekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe dekstre direkte al {destination}"
            },
            "slight left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}",
                "exit": "Direktiĝu al {exit} elveturejo ĉe maldekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe maldekstre direkte al {destination}"
            },
            "slight right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}",
                "exit": "Direktiĝu al {exit} elveturejo ĉe ldekstre",
                "exit_destination": "Direktiĝu al elveturejo {exit} ĉe dekstre direkte al {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Direktiĝu al enveturejo",
                "name": "Direktiĝu al enveturejo al {way_name}",
                "destination": "Direktiĝu al enveturejo direkte al {destination}"
            },
            "left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}"
            },
            "right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}"
            },
            "sharp left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}"
            },
            "sharp right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}"
            },
            "slight left": {
                "default": "Direktiĝu al enveturejo ĉe maldekstre",
                "name": "Direktiĝu al enveturejo ĉe maldekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe maldekstre al {destination}"
            },
            "slight right": {
                "default": "Direktiĝu al enveturejo ĉe dekstre",
                "name": "Direktiĝu al enveturejo ĉe dekstre al {way_name}",
                "destination": "Direktiĝu al enveturejo ĉe dekstre al {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Enveturu trafikcirklegon",
                    "name": "Enveturu trafikcirklegon kaj elveturu al {way_name}",
                    "destination": "Enveturu trafikcirklegon kaj elveturu direkte al {destination}"
                },
                "name": {
                    "default": "Enveturu {rotary_name}",
                    "name": "Enveturu {rotary_name} kaj elveturu al {way_name}",
                    "destination": "Enveturu {rotary_name} kaj elveturu direkte al {destination}"
                },
                "exit": {
                    "default": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo",
                    "name": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo al {way_name}",
                    "destination": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo direkte al {destination}"
                },
                "name_exit": {
                    "default": "Enveturu {rotary_name} kaj sekve al {exit_number} elveturejo",
                    "name": "Enveturu {rotary_name} kaj sekve al {exit_number} elveturejo al {way_name}",
                    "destination": "Enveturu {rotary_name} kaj sekve al {exit_number} elveturejo direkte al {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo",
                    "name": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo al {way_name}",
                    "destination": "Enveturu trafikcirklegon kaj sekve al {exit_number} elveturejo direkte al {destination}"
                },
                "default": {
                    "default": "Enveturu trafikcirklegon",
                    "name": "Enveturu trafikcirklegon kaj elveturu al {way_name}",
                    "destination": "Enveturu trafikcirklegon kaj elveturu direkte al {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Veturu {modifier}",
                "name": "Veturu {modifier} al {way_name}",
                "destination": "Veturu {modifier} direkte al {destination}"
            },
            "left": {
                "default": "Turniĝu maldekstren",
                "name": "Turniĝu maldekstren al {way_name}",
                "destination": "Turniĝu maldekstren direkte al {destination}"
            },
            "right": {
                "default": "Turniĝu dekstren",
                "name": "Turniĝu dekstren al {way_name}",
                "destination": "Turniĝu dekstren direkte al {destination}"
            },
            "straight": {
                "default": "Pluu rekten",
                "name": "Veturu rekten al {way_name}",
                "destination": "Veturu rekten direkte al {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Elveturu trafikcirklegon",
                "name": "Elveturu trafikcirklegon al {way_name}",
                "destination": "Elveturu trafikcirklegon direkte al {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Eliru trafikcirklegon",
                "name": "Elveturu trafikcirklegon al {way_name}",
                "destination": "Elveturu trafikcirklegon direkte al {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Veturu {modifier}",
                "name": "Veturu {modifier} al {way_name}",
                "destination": "Veturu {modifier} direkte al {destination}"
            },
            "left": {
                "default": "Turniĝu maldekstren",
                "name": "Turniĝu maldekstren al {way_name}",
                "destination": "Turniĝu maldekstren direkte al {destination}"
            },
            "right": {
                "default": "Turniĝu dekstren",
                "name": "Turniĝu dekstren al {way_name}",
                "destination": "Turniĝu dekstren direkte al {destination}"
            },
            "straight": {
                "default": "Veturu rekten",
                "name": "Veturu rekten al {way_name}",
                "destination": "Veturu rekten direkte al {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Pluu rekten"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],27:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1ª",
                "2": "2ª",
                "3": "3ª",
                "4": "4ª",
                "5": "5ª",
                "6": "6ª",
                "7": "7ª",
                "8": "8ª",
                "9": "9ª",
                "10": "10ª"
            },
            "direction": {
                "north": "norte",
                "northeast": "noreste",
                "east": "este",
                "southeast": "sureste",
                "south": "sur",
                "southwest": "suroeste",
                "west": "oeste",
                "northwest": "noroeste"
            },
            "modifier": {
                "left": "a la izquierda",
                "right": "a la derecha",
                "sharp left": "cerrada a la izquierda",
                "sharp right": "cerrada a la derecha",
                "slight left": "ligeramente a la izquierda",
                "slight right": "ligeramente a la derecha",
                "straight": "recto",
                "uturn": "cambio de sentido"
            },
            "lanes": {
                "xo": "Mantente a la derecha",
                "ox": "Mantente a la izquierda",
                "xox": "Mantente en el medio",
                "oxo": "Mantente a la izquierda o a la derecha"
            }
        },
        "modes": {
            "ferry": {
                "default": "Coge el ferry",
                "name": "Coge el ferry {way_name}",
                "destination": "Coge el ferry hacia {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} y luego en {distance}, {instruction_two}",
            "two linked": "{instruction_one} y luego {instruction_two}",
            "one in distance": "A {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "salida {exit}"
        },
        "arrive": {
            "default": {
                "default": "Has llegado a tu {nth} destino",
                "upcoming": "Vas a llegar a tu {nth} destino",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}"
            },
            "left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "sharp left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "sharp right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "slight right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "slight left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "straight": {
                "default": "Has llegado a tu {nth} destino, en frente",
                "upcoming": "Vas a llegar a tu {nth} destino, en frente",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, en frente"
            }
        },
        "continue": {
            "default": {
                "default": "Gire {modifier}",
                "name": "Cruce {modifier} en {way_name}",
                "destination": "Gire {modifier} hacia {destination}",
                "exit": "Gire {modifier} en {way_name}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa en {way_name}",
                "destination": "Continúa hacia {destination}",
                "distance": "Continúa recto por {distance}",
                "namedistance": "Continúa recto en {way_name} por {distance}"
            },
            "sharp left": {
                "default": "Gire a la izquierda",
                "name": "Gire a la izquierda en {way_name}",
                "destination": "Gire a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gire a la derecha",
                "name": "Gire a la derecha en {way_name}",
                "destination": "Gire a la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Gire a la izquierda",
                "name": "Doble levemente a la izquierda en {way_name}",
                "destination": "Gire a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Gire a la izquierda",
                "name": "Doble levemente a la derecha en {way_name}",
                "destination": "Gire a la izquierda hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido y continúa en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Dirígete al {direction}",
                "name": "Dirígete al {direction} por {way_name}",
                "namedistance": "Dirígete al {direction} en {way_name} por {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Al final de la calle gira {modifier}",
                "name": "Al final de la calle gira {modifier} por {way_name}",
                "destination": "Al final de la calle gira {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Al final de la calle continúa recto",
                "name": "Al final de la calle continúa recto por {way_name}",
                "destination": "Al final de la calle continúa recto hacia {destination}"
            },
            "uturn": {
                "default": "Al final de la calle haz un cambio de sentido",
                "name": "Al final de la calle haz un cambio de sentido en {way_name}",
                "destination": "Al final de la calle haz un cambio de sentido hacia {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "Mantente {modifier} en el cruce",
                "name": "Mantente {modifier} por {way_name}",
                "destination": "Mantente {modifier} hacia {destination}"
            },
            "slight left": {
                "default": "Mantente a la izquierda en el cruce",
                "name": "Mantente a la izquierda por {way_name}",
                "destination": "Mantente a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Mantente a la derecha en el cruce",
                "name": "Mantente a la derecha por {way_name}",
                "destination": "Mantente a la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Gira la izquierda en el cruce",
                "name": "Gira a la izquierda por {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gira a la derecha en el cruce",
                "name": "Gira a la derecha por {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Incorpórate {modifier}",
                "name": "Incorpórate {modifier} por {way_name}",
                "destination": "Incorpórate {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Incorpórate",
                "name": "Incorpórate por {way_name}",
                "destination": "Incorpórate hacia {destination}"
            },
            "slight left": {
                "default": "Incorpórate a la izquierda",
                "name": "Incorpórate a la izquierda por {way_name}",
                "destination": "Incorpórate a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Incorpórate a la derecha",
                "name": "Incorpórate a la derecha por {way_name}",
                "destination": "Incorpórate a la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Incorpórate a la izquierda",
                "name": "Incorpórate a la izquierda por {way_name}",
                "destination": "Incorpórate a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Incorpórate a la derecha",
                "name": "Incorpórate a la derecha por {way_name}",
                "destination": "Incorpórate a la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continúa {modifier}",
                "name": "Continúa {modifier} por {way_name}",
                "destination": "Continúa {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa por {way_name}",
                "destination": "Continúa hacia {destination}"
            },
            "sharp left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda por {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha por {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Continúa ligeramente por la izquierda",
                "name": "Continúa ligeramente por la izquierda por {way_name}",
                "destination": "Continúa ligeramente por la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Continúa ligeramente por la derecha",
                "name": "Continúa ligeramente por la derecha por {way_name}",
                "destination": "Continúa ligeramente por la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continúa {modifier}",
                "name": "Continúa {modifier} por {way_name}",
                "destination": "Continúa {modifier} hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Coge la cuesta abajo",
                "name": "Coge la cuesta abajo por {way_name}",
                "destination": "Coge la cuesta abajo hacia {destination}",
                "exit": "Coge la cuesta abajo {exit}",
                "exit_destination": "Coge la cuesta abajo {exit} hacia {destination}"
            },
            "left": {
                "default": "Coge la cuesta abajo de la izquierda",
                "name": "Coge la cuesta abajo de la izquierda por {way_name}",
                "destination": "Coge la cuesta abajo de la izquierda hacia {destination}",
                "exit": "Coge la cuesta abajo {exit} a tu izquierda",
                "exit_destination": "Coge la cuesta abajo {exit} a tu izquierda hacia {destination}"
            },
            "right": {
                "default": "Coge la cuesta abajo de la derecha",
                "name": "Coge la cuesta abajo de la derecha por {way_name}",
                "destination": "Coge la cuesta abajo de la derecha hacia {destination}",
                "exit": "Coge la cuesta abajo {exit}",
                "exit_destination": "Coge la cuesta abajo {exit} hacia {destination}"
            },
            "sharp left": {
                "default": "Coge la cuesta abajo de la izquierda",
                "name": "Coge la cuesta abajo de la izquierda por {way_name}",
                "destination": "Coge la cuesta abajo de la izquierda hacia {destination}",
                "exit": "Coge la cuesta abajo {exit} a tu izquierda",
                "exit_destination": "Coge la cuesta abajo {exit} a tu izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Coge la cuesta abajo de la derecha",
                "name": "Coge la cuesta abajo de la derecha por {way_name}",
                "destination": "Coge la cuesta abajo de la derecha hacia {destination}",
                "exit": "Coge la cuesta abajo {exit}",
                "exit_destination": "Coge la cuesta abajo {exit} hacia {destination}"
            },
            "slight left": {
                "default": "Coge la cuesta abajo de la izquierda",
                "name": "Coge la cuesta abajo de la izquierda por {way_name}",
                "destination": "Coge la cuesta abajo de la izquierda hacia {destination}",
                "exit": "Coge la cuesta abajo {exit} a tu izquierda",
                "exit_destination": "Coge la cuesta abajo {exit} a tu izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Coge la cuesta abajo de la derecha",
                "name": "Coge la cuesta abajo de la derecha por {way_name}",
                "destination": "Coge la cuesta abajo de la derecha hacia {destination}",
                "exit": "Coge la cuesta abajo {exit}",
                "exit_destination": "Coge la cuesta abajo {exit} hacia {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Coge la cuesta",
                "name": "Coge la cuesta por {way_name}",
                "destination": "Coge la cuesta hacia {destination}"
            },
            "left": {
                "default": "Coge la cuesta de la izquierda",
                "name": "Coge la cuesta de la izquierda por {way_name}",
                "destination": "Coge la cuesta de la izquierda hacia {destination}"
            },
            "right": {
                "default": "Coge la cuesta de la derecha",
                "name": "Coge la cuesta de la derecha por {way_name}",
                "destination": "Coge la cuesta de la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Coge la cuesta de la izquierda",
                "name": "Coge la cuesta de la izquierda por {way_name}",
                "destination": "Coge la cuesta de la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Coge la cuesta de la derecha",
                "name": "Coge la cuesta de la derecha por {way_name}",
                "destination": "Coge la cuesta de la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Coge la cuesta de la izquierda",
                "name": "Coge la cuesta de la izquierda por {way_name}",
                "destination": "Coge la cuesta de la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Coge la cuesta de la derecha",
                "name": "Coge la cuesta de la derecha por {way_name}",
                "destination": "Coge la cuesta de la derecha hacia {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Incorpórate en la rotonda",
                    "name": "En la rotonda sal por {way_name}",
                    "destination": "En la rotonda sal hacia {destination}"
                },
                "name": {
                    "default": "En {rotary_name}",
                    "name": "En {rotary_name} sal por {way_name}",
                    "destination": "En {rotary_name} sal hacia {destination}"
                },
                "exit": {
                    "default": "En la rotonda toma la {exit_number} salida",
                    "name": "En la rotonda toma la {exit_number} salida por {way_name}",
                    "destination": "En la rotonda toma la {exit_number} salida hacia {destination}"
                },
                "name_exit": {
                    "default": "En {rotary_name} toma la {exit_number} salida",
                    "name": "En {rotary_name} toma la {exit_number} salida por {way_name}",
                    "destination": "En {rotary_name} toma la {exit_number} salida hacia {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "En la rotonda toma la {exit_number} salida",
                    "name": "En la rotonda toma la {exit_number} salida por {way_name}",
                    "destination": "En la rotonda toma la {exit_number} salida hacia {destination}"
                },
                "default": {
                    "default": "Incorpórate en la rotonda",
                    "name": "Incorpórate en la rotonda y sal en {way_name}",
                    "destination": "Incorpórate en la rotonda y sal hacia {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Siga {modifier}",
                "name": "Siga {modifier} en {way_name}",
                "destination": "Siga {modifier} hacia {destination}"
            },
            "left": {
                "default": "Gire a la izquierda",
                "name": "Gire a la izquierda en {way_name}",
                "destination": "Gire a la izquierda hacia {destination}"
            },
            "right": {
                "default": "Gire a la derecha",
                "name": "Gire a la derecha en {way_name}",
                "destination": "Gire a la derecha hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa recto por {way_name}",
                "destination": "Continúa recto hacia {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Sal la rotonda",
                "name": "Toma la salida por {way_name}",
                "destination": "Toma la salida hacia {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Sal la rotonda",
                "name": "Toma la salida por {way_name}",
                "destination": "Toma la salida hacia {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Gira {modifier}",
                "name": "Gira {modifier} por {way_name}",
                "destination": "Gira {modifier} hacia {destination}"
            },
            "left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda por {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha por {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa recto por {way_name}",
                "destination": "Continúa recto hacia {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continúa recto"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],28:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1ª",
                "2": "2ª",
                "3": "3ª",
                "4": "4ª",
                "5": "5ª",
                "6": "6ª",
                "7": "7ª",
                "8": "8ª",
                "9": "9ª",
                "10": "10ª"
            },
            "direction": {
                "north": "norte",
                "northeast": "noreste",
                "east": "este",
                "southeast": "sureste",
                "south": "sur",
                "southwest": "suroeste",
                "west": "oeste",
                "northwest": "noroeste"
            },
            "modifier": {
                "left": "izquierda",
                "right": "derecha",
                "sharp left": "cerrada a la izquierda",
                "sharp right": "cerrada a la derecha",
                "slight left": "levemente a la izquierda",
                "slight right": "levemente a la derecha",
                "straight": "recto",
                "uturn": "cambio de sentido"
            },
            "lanes": {
                "xo": "Mantente a la derecha",
                "ox": "Mantente a la izquierda",
                "xox": "Mantente en el medio",
                "oxo": "Mantente a la izquierda o derecha"
            }
        },
        "modes": {
            "ferry": {
                "default": "Coge el ferry",
                "name": "Coge el ferry {way_name}",
                "destination": "Coge el ferry a {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} y luego a {distance}, {instruction_two}",
            "two linked": "{instruction_one} y luego {instruction_two}",
            "one in distance": "A {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "salida {exit}"
        },
        "arrive": {
            "default": {
                "default": "Has llegado a tu {nth} destino",
                "upcoming": "Vas a llegar a tu {nth} destino",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}"
            },
            "left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "sharp left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "sharp right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "slight right": {
                "default": "Has llegado a tu {nth} destino, a la derecha",
                "upcoming": "Vas a llegar a tu {nth} destino, a la derecha",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la derecha"
            },
            "slight left": {
                "default": "Has llegado a tu {nth} destino, a la izquierda",
                "upcoming": "Vas a llegar a tu {nth} destino, a la izquierda",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, a la izquierda"
            },
            "straight": {
                "default": "Has llegado a tu {nth} destino, en frente",
                "upcoming": "Vas a llegar a tu {nth} destino, en frente",
                "short": "Has llegado",
                "short-upcoming": "Vas a llegar",
                "named": "Has llegado a {waypoint_name}, en frente"
            }
        },
        "continue": {
            "default": {
                "default": "Gira a {modifier}",
                "name": "Cruza a la{modifier}  en {way_name}",
                "destination": "Gira a {modifier} hacia {destination}",
                "exit": "Gira a {modifier} en {way_name}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa en {way_name}",
                "destination": "Continúa hacia {destination}",
                "distance": "Continúa recto por {distance}",
                "namedistance": "Continúa recto en {way_name} por {distance}"
            },
            "sharp left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha en {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Gira a la izquierda",
                "name": "Dobla levemente a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Gira a la izquierda",
                "name": "Dobla levemente a la derecha en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido y continúa en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Ve a {direction}",
                "name": "Ve a {direction} en {way_name}",
                "namedistance": "Ve a {direction} en {way_name} por {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Gira  a {modifier}",
                "name": "Gira a {modifier} en {way_name}",
                "destination": "Gira a {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa recto en {way_name}",
                "destination": "Continúa recto hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido al final de la via",
                "name": "Haz un cambio de sentido en {way_name} al final de la via",
                "destination": "Haz un cambio de sentido hacia {destination} al final de la via"
            }
        },
        "fork": {
            "default": {
                "default": "Mantente  {modifier} en el cruza",
                "name": "Mantente {modifier} en {way_name}",
                "destination": "Mantente {modifier} hacia {destination}"
            },
            "slight left": {
                "default": "Mantente a la izquierda en el cruza",
                "name": "Mantente a la izquierda en {way_name}",
                "destination": "Mantente a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Mantente a la derecha en el cruza",
                "name": "Mantente a la derecha en {way_name}",
                "destination": "Mantente a la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Gira a la izquierda en el cruza",
                "name": "Gira a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gira a la derecha en el cruza",
                "name": "Gira a la derecha en {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Incorpórate a {modifier}",
                "name": "Incorpórate a {modifier} en {way_name}",
                "destination": "Incorpórate a {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Incorpórate",
                "name": "Incorpórate a {way_name}",
                "destination": "Incorpórate hacia {destination}"
            },
            "slight left": {
                "default": "Incorpórate a la izquierda",
                "name": "Incorpórate a la izquierda en {way_name}",
                "destination": "Incorpórate a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Incorpórate a la derecha",
                "name": "Incorpórate a la derecha en {way_name}",
                "destination": "Incorpórate a la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Incorpórate a la izquierda",
                "name": "Incorpórate a la izquierda en {way_name}",
                "destination": "Incorpórate a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Incorpórate a la derecha",
                "name": "Incorpórate a la derecha en {way_name}",
                "destination": "Incorpórate a la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continúa {modifier}",
                "name": "Continúa {modifier} en {way_name}",
                "destination": "Continúa {modifier} hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa en {way_name}",
                "destination": "Continúa hacia {destination}"
            },
            "sharp left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha en {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Continúa levemente a la izquierda",
                "name": "Continúa levemente a la izquierda en {way_name}",
                "destination": "Continúa levemente a la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Continúa levemente a la derecha",
                "name": "Continúa levemente a la derecha en {way_name}",
                "destination": "Continúa levemente a la derecha hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continúa {modifier}",
                "name": "Continúa {modifier} en {way_name}",
                "destination": "Continúa {modifier} hacia {destination}"
            },
            "uturn": {
                "default": "Haz un cambio de sentido",
                "name": "Haz un cambio de sentido en {way_name}",
                "destination": "Haz un cambio de sentido hacia {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Toma la salida",
                "name": "Toma la salida en {way_name}",
                "destination": "Toma la salida hacia {destination}",
                "exit": "Toma la salida {exit}",
                "exit_destination": "Toma la salida {exit} hacia {destination}"
            },
            "left": {
                "default": "Toma la salida en la izquierda",
                "name": "Toma la salida en la izquierda en {way_name}",
                "destination": "Toma la salida en la izquierda en {destination}",
                "exit": "Toma la salida {exit} en la izquierda",
                "exit_destination": "Toma la salida {exit} en la izquierda hacia {destination}"
            },
            "right": {
                "default": "Toma la salida en la derecha",
                "name": "Toma la salida en la derecha en {way_name}",
                "destination": "Toma la salida en la derecha hacia {destination}",
                "exit": "Toma la salida {exit} en la derecha",
                "exit_destination": "Toma la salida {exit} en la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Ve cuesta abajo en la izquierda",
                "name": "Ve cuesta abajo en la izquierda en {way_name}",
                "destination": "Ve cuesta abajo en la izquierda hacia {destination}",
                "exit": "Toma la salida {exit} en la izquierda",
                "exit_destination": "Toma la salida {exit} en la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Ve cuesta abajo en la derecha",
                "name": "Ve cuesta abajo en la derecha en {way_name}",
                "destination": "Ve cuesta abajo en la derecha hacia {destination}",
                "exit": "Toma la salida {exit} en la derecha",
                "exit_destination": "Toma la salida {exit} en la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Ve cuesta abajo en la izquierda",
                "name": "Ve cuesta abajo en la izquierda en {way_name}",
                "destination": "Ve cuesta abajo en la izquierda hacia {destination}",
                "exit": "Toma la salida {exit} en la izquierda",
                "exit_destination": "Toma la salida {exit} en la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Toma la salida en la derecha",
                "name": "Toma la salida en la derecha en {way_name}",
                "destination": "Toma la salida en la derecha hacia {destination}",
                "exit": "Toma la salida {exit} en la derecha",
                "exit_destination": "Toma la salida {exit} en la derecha hacia {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Toma la rampa",
                "name": "Toma la rampa en {way_name}",
                "destination": "Toma la rampa hacia {destination}"
            },
            "left": {
                "default": "Toma la rampa en la izquierda",
                "name": "Toma la rampa en la izquierda en {way_name}",
                "destination": "Toma la rampa en la izquierda hacia {destination}"
            },
            "right": {
                "default": "Toma la rampa en la derecha",
                "name": "Toma la rampa en la derecha en {way_name}",
                "destination": "Toma la rampa en la derecha hacia {destination}"
            },
            "sharp left": {
                "default": "Toma la rampa en la izquierda",
                "name": "Toma la rampa en la izquierda en {way_name}",
                "destination": "Toma la rampa en la izquierda hacia {destination}"
            },
            "sharp right": {
                "default": "Toma la rampa en la derecha",
                "name": "Toma la rampa en la derecha en {way_name}",
                "destination": "Toma la rampa en la derecha hacia {destination}"
            },
            "slight left": {
                "default": "Toma la rampa en la izquierda",
                "name": "Toma la rampa en la izquierda en {way_name}",
                "destination": "Toma la rampa en la izquierda hacia {destination}"
            },
            "slight right": {
                "default": "Toma la rampa en la derecha",
                "name": "Toma la rampa en la derecha en {way_name}",
                "destination": "Toma la rampa en la derecha hacia {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Entra en la rotonda",
                    "name": "Entra en la rotonda y sal en {way_name}",
                    "destination": "Entra en la rotonda y sal hacia {destination}"
                },
                "name": {
                    "default": "Entra en {rotary_name}",
                    "name": "Entra en {rotary_name} y sal en {way_name}",
                    "destination": "Entra en {rotary_name} y sal hacia {destination}"
                },
                "exit": {
                    "default": "Entra en la rotonda y toma la {exit_number} salida",
                    "name": "Entra en la rotonda y toma la {exit_number} salida a {way_name}",
                    "destination": "Entra en la rotonda y toma la {exit_number} salida hacia {destination}"
                },
                "name_exit": {
                    "default": "Entra en {rotary_name} y coge la {exit_number} salida",
                    "name": "Entra en {rotary_name} y coge la {exit_number} salida en {way_name}",
                    "destination": "Entra en {rotary_name} y coge la {exit_number} salida hacia {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Entra en la rotonda y toma la {exit_number} salida",
                    "name": "Entra en la rotonda y toma la {exit_number} salida a {way_name}",
                    "destination": "Entra en la rotonda y toma la {exit_number} salida hacia {destination}"
                },
                "default": {
                    "default": "Entra en la rotonda",
                    "name": "Entra en la rotonda y sal en {way_name}",
                    "destination": "Entra en la rotonda y sal hacia {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Sigue {modifier}",
                "name": "Sigue {modifier} en {way_name}",
                "destination": "Sigue {modifier} hacia {destination}"
            },
            "left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha en {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "straight": {
                "default": "Continúa recto",
                "name": "Continúa recto en {way_name}",
                "destination": "Continúa recto hacia {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Sal la rotonda",
                "name": "Sal la rotonda en {way_name}",
                "destination": "Sal la rotonda hacia {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Sal la rotonda",
                "name": "Sal la rotonda en {way_name}",
                "destination": "Sal la rotonda hacia {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Sigue {modifier}",
                "name": "Sigue {modifier} en {way_name}",
                "destination": "Sigue {modifier} hacia {destination}"
            },
            "left": {
                "default": "Gira a la izquierda",
                "name": "Gira a la izquierda en {way_name}",
                "destination": "Gira a la izquierda hacia {destination}"
            },
            "right": {
                "default": "Gira a la derecha",
                "name": "Gira a la derecha en {way_name}",
                "destination": "Gira a la derecha hacia {destination}"
            },
            "straight": {
                "default": "Ve recto",
                "name": "Ve recto en {way_name}",
                "destination": "Ve recto hacia {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continúa recto"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],29:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1.",
                "2": "2.",
                "3": "3.",
                "4": "4.",
                "5": "5.",
                "6": "6.",
                "7": "7.",
                "8": "8.",
                "9": "9.",
                "10": "10."
            },
            "direction": {
                "north": "pohjoiseen",
                "northeast": "koilliseen",
                "east": "itään",
                "southeast": "kaakkoon",
                "south": "etelään",
                "southwest": "lounaaseen",
                "west": "länteen",
                "northwest": "luoteeseen"
            },
            "modifier": {
                "left": "vasemmall(e/a)",
                "right": "oikeall(e/a)",
                "sharp left": "jyrkästi vasempaan",
                "sharp right": "jyrkästi oikeaan",
                "slight left": "loivasti vasempaan",
                "slight right": "loivasti oikeaan",
                "straight": "suoraan eteenpäin",
                "uturn": "U-käännös"
            },
            "lanes": {
                "xo": "Pysy oikealla",
                "ox": "Pysy vasemmalla",
                "xox": "Pysy keskellä",
                "oxo": "Pysy vasemmalla tai oikealla"
            }
        },
        "modes": {
            "ferry": {
                "default": "Aja lautalle",
                "name": "Aja lautalle {way_name}",
                "destination": "Aja lautalle, jonka määränpää on {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, sitten {distance} päästä, {instruction_two}",
            "two linked": "{instruction_one}, sitten {instruction_two}",
            "one in distance": "{distance} päästä, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "{exit}"
        },
        "arrive": {
            "default": {
                "default": "Olet saapunut {nth} määränpäähäsi",
                "upcoming": "Saavut {nth} määränpäähäsi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}"
            },
            "left": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on vasemmalla puolellasi"
            },
            "right": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on oikealla puolellasi"
            },
            "sharp left": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on vasemmalla puolellasi"
            },
            "sharp right": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on oikealla puolellasi"
            },
            "slight right": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on oikealla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on oikealla puolellasi"
            },
            "slight left": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "upcoming": "Saavut {nth} määränpäähäsi, joka on vasemmalla puolellasi",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on vasemmalla puolellasi"
            },
            "straight": {
                "default": "Olet saapunut {nth} määränpäähäsi, joka on suoraan edessäsi",
                "upcoming": "Saavut {nth} määränpäähäsi, suoraan edessä",
                "short": "Olet saapunut",
                "short-upcoming": "Saavut",
                "named": "Olet saapunut määränpäähän {waypoint_name}, joka on suoraan edessäsi"
            }
        },
        "continue": {
            "default": {
                "default": "Käänny {modifier}",
                "name": "Käänny {modifier} pysyäksesi tiellä {way_name}",
                "destination": "Käänny {modifier} suuntana {destination}",
                "exit": "Käänny {modifier} tielle {way_name}"
            },
            "straight": {
                "default": "Jatka suoraan eteenpäin",
                "name": "Jatka suoraan pysyäksesi tiellä {way_name}",
                "destination": "Jatka suuntana {destination}",
                "distance": "Jatka suoraan {distance}",
                "namedistance": "Jatka tiellä {way_name} {distance}"
            },
            "sharp left": {
                "default": "Jatka jyrkästi vasempaan",
                "name": "Jatka jyrkästi vasempaan pysyäksesi tiellä {way_name}",
                "destination": "Jatka jyrkästi vasempaan suuntana {destination}"
            },
            "sharp right": {
                "default": "Jatka jyrkästi oikeaan",
                "name": "Jatka jyrkästi oikeaan pysyäksesi tiellä {way_name}",
                "destination": "Jatka jyrkästi oikeaan suuntana {destination}"
            },
            "slight left": {
                "default": "Jatka loivasti vasempaan",
                "name": "Jatka loivasti vasempaan pysyäksesi tiellä {way_name}",
                "destination": "Jatka loivasti vasempaan suuntana {destination}"
            },
            "slight right": {
                "default": "Jatka loivasti oikeaan",
                "name": "Jatka loivasti oikeaan pysyäksesi tiellä {way_name}",
                "destination": "Jatka loivasti oikeaan suuntana {destination}"
            },
            "uturn": {
                "default": "Tee U-käännös",
                "name": "Tee U-käännös ja jatka tietä {way_name}",
                "destination": "Tee U-käännös suuntana {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Aja {direction}",
                "name": "Aja tietä {way_name} {direction}",
                "namedistance": "Aja {distance} {direction} tietä {way_name} "
            }
        },
        "end of road": {
            "default": {
                "default": "Käänny {modifier}",
                "name": "Käänny {modifier} tielle {way_name}",
                "destination": "Käänny {modifier} suuntana {destination}"
            },
            "straight": {
                "default": "Jatka suoraan eteenpäin",
                "name": "Jatka suoraan eteenpäin tielle {way_name}",
                "destination": "Jatka suoraan eteenpäin suuntana {destination}"
            },
            "uturn": {
                "default": "Tien päässä tee U-käännös",
                "name": "Tien päässä tee U-käännös tielle {way_name}",
                "destination": "Tien päässä tee U-käännös suuntana {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "Jatka tienhaarassa {modifier}",
                "name": "Jatka {modifier} tielle {way_name}",
                "destination": "Jatka {modifier} suuntana {destination}"
            },
            "slight left": {
                "default": "Pysy vasemmalla tienhaarassa",
                "name": "Pysy vasemmalla tielle {way_name}",
                "destination": "Pysy vasemmalla suuntana {destination}"
            },
            "slight right": {
                "default": "Pysy oikealla tienhaarassa",
                "name": "Pysy oikealla tielle {way_name}",
                "destination": "Pysy oikealla suuntana {destination}"
            },
            "sharp left": {
                "default": "Käänny tienhaarassa jyrkästi vasempaan",
                "name": "Käänny tienhaarassa jyrkästi vasempaan tielle {way_name}",
                "destination": "Käänny tienhaarassa jyrkästi vasempaan suuntana {destination}"
            },
            "sharp right": {
                "default": "Käänny tienhaarassa jyrkästi oikeaan",
                "name": "Käänny tienhaarassa jyrkästi oikeaan tielle {way_name}",
                "destination": "Käänny tienhaarassa jyrkästi oikeaan suuntana {destination}"
            },
            "uturn": {
                "default": "Tee U-käännös",
                "name": "Tee U-käännös tielle {way_name}",
                "destination": "Tee U-käännös suuntana {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Liity {modifier}",
                "name": "Liity {modifier}, tielle {way_name}",
                "destination": "Liity {modifier}, suuntana {destination}"
            },
            "straight": {
                "default": "Liity",
                "name": "Liity tielle {way_name}",
                "destination": "Liity suuntana {destination}"
            },
            "slight left": {
                "default": "Liity vasemmalle",
                "name": "Liity vasemmalle, tielle {way_name}",
                "destination": "Liity vasemmalle, suuntana {destination}"
            },
            "slight right": {
                "default": "Liity oikealle",
                "name": "Liity oikealle, tielle {way_name}",
                "destination": "Liity oikealle, suuntana {destination}"
            },
            "sharp left": {
                "default": "Liity vasemmalle",
                "name": "Liity vasemmalle, tielle {way_name}",
                "destination": "Liity vasemmalle, suuntana {destination}"
            },
            "sharp right": {
                "default": "Liity oikealle",
                "name": "Liity oikealle, tielle {way_name}",
                "destination": "Liity oikealle, suuntana {destination}"
            },
            "uturn": {
                "default": "Tee U-käännös",
                "name": "Tee U-käännös tielle {way_name}",
                "destination": "Tee U-käännös suuntana {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Jatka {modifier}",
                "name": "Jatka {modifier} tielle {way_name}",
                "destination": "Jatka {modifier} suuntana {destination}"
            },
            "straight": {
                "default": "Jatka suoraan eteenpäin",
                "name": "Jatka tielle {way_name}",
                "destination": "Jatka suuntana {destination}"
            },
            "sharp left": {
                "default": "Käänny jyrkästi vasempaan",
                "name": "Käänny jyrkästi vasempaan tielle {way_name}",
                "destination": "Käänny jyrkästi vasempaan suuntana {destination}"
            },
            "sharp right": {
                "default": "Käänny jyrkästi oikeaan",
                "name": "Käänny jyrkästi oikeaan tielle {way_name}",
                "destination": "Käänny jyrkästi oikeaan suuntana {destination}"
            },
            "slight left": {
                "default": "Jatka loivasti vasempaan",
                "name": "Jatka loivasti vasempaan tielle {way_name}",
                "destination": "Jatka loivasti vasempaan suuntana {destination}"
            },
            "slight right": {
                "default": "Jatka loivasti oikeaan",
                "name": "Jatka loivasti oikeaan tielle {way_name}",
                "destination": "Jatka loivasti oikeaan suuntana {destination}"
            },
            "uturn": {
                "default": "Tee U-käännös",
                "name": "Tee U-käännös tielle {way_name}",
                "destination": "Tee U-käännös suuntana {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Jatka {modifier}",
                "name": "Jatka {modifier} tielle {way_name}",
                "destination": "Jatka {modifier} suuntana {destination}"
            },
            "uturn": {
                "default": "Tee U-käännös",
                "name": "Tee U-käännös tielle {way_name}",
                "destination": "Tee U-käännös suuntana {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Aja erkanemiskaistalle",
                "name": "Aja erkanemiskaistaa tielle {way_name}",
                "destination": "Aja erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit}",
                "exit_destination": "Ota poistuminen {exit}, suuntana {destination}"
            },
            "left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} vasemmalla",
                "exit_destination": "Ota poistuminen {exit} vasemmalla, suuntana {destination}"
            },
            "right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} oikealla",
                "exit_destination": "Ota poistuminen {exit} oikealla, suuntana {destination}"
            },
            "sharp left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} vasemmalla",
                "exit_destination": "Ota poistuminen {exit} vasemmalla, suuntana {destination}"
            },
            "sharp right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} oikealla",
                "exit_destination": "Ota poistuminen {exit} oikealla, suuntana {destination}"
            },
            "slight left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} vasemmalla",
                "exit_destination": "Ota poistuminen {exit} vasemmalla, suuntana {destination}"
            },
            "slight right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}",
                "exit": "Ota poistuminen {exit} oikealla",
                "exit_destination": "Ota poistuminen {exit} oikealla, suuntana {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Aja erkanemiskaistalle",
                "name": "Aja erkanemiskaistaa tielle {way_name}",
                "destination": "Aja erkanemiskaistalle suuntana {destination}"
            },
            "left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}"
            },
            "right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}"
            },
            "sharp left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}"
            },
            "sharp right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}"
            },
            "slight left": {
                "default": "Aja vasemmalla olevalle erkanemiskaistalle",
                "name": "Aja vasemmalla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja vasemmalla olevalle erkanemiskaistalle suuntana {destination}"
            },
            "slight right": {
                "default": "Aja oikealla olevalle erkanemiskaistalle",
                "name": "Aja oikealla olevaa erkanemiskaistaa tielle {way_name}",
                "destination": "Aja oikealla olevalle erkanemiskaistalle suuntana {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Aja liikenneympyrään",
                    "name": "Aja liikenneympyrään ja valitse erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään ja valitse erkanemiskaista suuntana {destination}"
                },
                "name": {
                    "default": "Aja liikenneympyrään {rotary_name}",
                    "name": "Aja liikenneympyrään {rotary_name} ja valitse erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään {rotary_name} ja valitse erkanemiskaista suuntana {destination}"
                },
                "exit": {
                    "default": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista",
                    "name": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista suuntana {destination}"
                },
                "name_exit": {
                    "default": "Aja liikenneympyrään {rotary_name} ja valitse {exit_number} erkanemiskaista",
                    "name": "Aja liikenneympyrään {rotary_name} ja valitse {exit_number} erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään {rotary_name} ja valitse {exit_number} erkanemiskaista suuntana {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista",
                    "name": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään ja valitse {exit_number} erkanemiskaista suuntana {destination}"
                },
                "default": {
                    "default": "Aja liikenneympyrään",
                    "name": "Aja liikenneympyrään ja valitse erkanemiskaista tielle {way_name}",
                    "destination": "Aja liikenneympyrään ja valitse erkanemiskaista suuntana {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Käänny {modifier}",
                "name": "Käänny {modifier} tielle {way_name}",
                "destination": "Käänny {modifier} suuntana {destination}"
            },
            "left": {
                "default": "Käänny vasempaan",
                "name": "Käänny vasempaan tielle {way_name}",
                "destination": "Käänny vasempaan suuntana {destination}"
            },
            "right": {
                "default": "Käänny oikeaan",
                "name": "Käänny oikeaan tielle {way_name}",
                "destination": "Käänny oikeaan suuntana {destination}"
            },
            "straight": {
                "default": "Jatka suoraan eteenpäin",
                "name": "Jatka suoraan eteenpäin tielle {way_name}",
                "destination": "Jatka suoraan eteenpäin suuntana {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Poistu liikenneympyrästä",
                "name": "Poistu liikenneympyrästä tielle {way_name}",
                "destination": "Poistu liikenneympyrästä suuntana {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Poistu liikenneympyrästä",
                "name": "Poistu liikenneympyrästä tielle {way_name}",
                "destination": "Poistu liikenneympyrästä suuntana {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Käänny {modifier}",
                "name": "Käänny {modifier} tielle {way_name}",
                "destination": "Käänny {modifier} suuntana {destination}"
            },
            "left": {
                "default": "Käänny vasempaan",
                "name": "Käänny vasempaan tielle {way_name}",
                "destination": "Käänny vasempaan suuntana {destination}"
            },
            "right": {
                "default": "Käänny oikeaan",
                "name": "Käänny oikeaan tielle {way_name}",
                "destination": "Käänny oikeaan suuntana {destination}"
            },
            "straight": {
                "default": "Aja suoraan eteenpäin",
                "name": "Aja suoraan eteenpäin tielle {way_name}",
                "destination": "Aja suoraan eteenpäin suuntana {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Jatka suoraan eteenpäin"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],30:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "première",
                "2": "seconde",
                "3": "troisième",
                "4": "quatrième",
                "5": "cinquième",
                "6": "sixième",
                "7": "septième",
                "8": "huitième",
                "9": "neuvième",
                "10": "dixième"
            },
            "direction": {
                "north": "le nord",
                "northeast": "le nord-est",
                "east": "l’est",
                "southeast": "le sud-est",
                "south": "le sud",
                "southwest": "le sud-ouest",
                "west": "l’ouest",
                "northwest": "le nord-ouest"
            },
            "modifier": {
                "left": "à gauche",
                "right": "à droite",
                "sharp left": "franchement à gauche",
                "sharp right": "franchement à droite",
                "slight left": "légèrement à gauche",
                "slight right": "légèrement à droite",
                "straight": "tout droit",
                "uturn": "demi-tour"
            },
            "lanes": {
                "xo": "Tenir la droite",
                "ox": "Tenir la gauche",
                "xox": "Rester au milieu",
                "oxo": "Tenir la gauche ou la droite"
            }
        },
        "modes": {
            "ferry": {
                "default": "Prendre le ferry",
                "name": "Prendre le ferry {way_name:article}",
                "destination": "Prendre le ferry en direction {destination:preposition}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, puis, dans {distance}, {instruction_two}",
            "two linked": "{instruction_one}, puis {instruction_two}",
            "one in distance": "Dans {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "sortie n°{exit}"
        },
        "arrive": {
            "default": {
                "default": "Vous êtes arrivé à votre {nth} destination",
                "upcoming": "Vous arriverez à votre {nth} destination",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé {waypoint_name:arrival}"
            },
            "left": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la gauche",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la gauche",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, sur la gauche"
            },
            "right": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la droite",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la droite",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé à  {waypoint_name:arrival}, sur la droite"
            },
            "sharp left": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la gauche",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la gauche",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, sur la gauche"
            },
            "sharp right": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la droite",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la droite",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, sur la droite"
            },
            "slight right": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la droite",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la droite",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous arriverez",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, sur la droite"
            },
            "slight left": {
                "default": "Vous êtes arrivé à votre {nth} destination, sur la gauche",
                "upcoming": "Vous arriverez à votre {nth} destination, sur la gauche",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous êtes arrivé",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, sur la gauche"
            },
            "straight": {
                "default": "Vous êtes arrivé à votre {nth} destination, droit devant",
                "upcoming": "Vous arriverez à votre {nth} destination, droit devant",
                "short": "Vous êtes arrivé",
                "short-upcoming": "Vous êtes arrivé",
                "named": "Vous êtes arrivé {waypoint_name:arrival}, droit devant"
            }
        },
        "continue": {
            "default": {
                "default": "Tourner {modifier}",
                "name": "Tourner {modifier} pour rester sur {way_name:article}",
                "destination": "Tourner {modifier} en direction {destination:preposition}",
                "exit": "Tourner {modifier} sur {way_name:article}"
            },
            "straight": {
                "default": "Continuer tout droit",
                "name": "Continuer tout droit pour rester sur {way_name:article}",
                "destination": "Continuer tout droit en direction {destination:preposition}",
                "distance": "Continuer tout droit sur {distance}",
                "namedistance": "Continuer sur {way_name:article} sur {distance}"
            },
            "sharp left": {
                "default": "Tourner franchement à gauche",
                "name": "Tourner franchement à gauche pour rester sur {way_name:article}",
                "destination": "Tourner franchement à gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "Tourner franchement à droite",
                "name": "Tourner franchement à droite pour rester sur {way_name:article}",
                "destination": "Tourner franchement à droite en direction {destination:preposition}"
            },
            "slight left": {
                "default": "Tourner légèrement à gauche",
                "name": "Tourner légèrement à gauche pour rester sur {way_name:article}",
                "destination": "Tourner légèrement à gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "Tourner légèrement à droite",
                "name": "Tourner légèrement à droite pour rester sur {way_name:article}",
                "destination": "Tourner légèrement à droite en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour",
                "name": "Faire demi-tour et continuer sur {way_name:article}",
                "destination": "Faire demi-tour en direction {destination:preposition}"
            }
        },
        "depart": {
            "default": {
                "default": "Se diriger vers {direction}",
                "name": "Se diriger vers {direction} sur {way_name:article}",
                "namedistance": "Se diriger vers {direction} sur {way_name:article} sur {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Tourner {modifier}",
                "name": "Tourner {modifier} sur {way_name:article}",
                "destination": "Tourner {modifier} en direction {destination:preposition}"
            },
            "straight": {
                "default": "Continuer tout droit",
                "name": "Continuer tout droit sur {way_name:article}",
                "destination": "Continuer tout droit en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour à la fin de la route",
                "name": "Faire demi-tour à la fin {way_name:preposition}",
                "destination": "Faire demi-tour à la fin de la route en direction {destination:preposition}"
            }
        },
        "fork": {
            "default": {
                "default": "Tenir {modifier} à l’embranchement",
                "name": "Tenir {modifier} sur {way_name:article}",
                "destination": "Tenir {modifier} en direction {destination:preposition}"
            },
            "slight left": {
                "default": "Tenir la gauche à l’embranchement",
                "name": "Tenir la gauche sur {way_name:article}",
                "destination": "Tenir la gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "Tenir la droite à l’embranchement",
                "name": "Tenir la droite sur {way_name:article}",
                "destination": "Tenir la droite en direction {destination:preposition}"
            },
            "sharp left": {
                "default": "Tourner franchement à gauche à l’embranchement",
                "name": "Tourner franchement à gauche sur {way_name:article}",
                "destination": "Tourner franchement à gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "Tourner franchement à droite à l’embranchement",
                "name": "Tourner franchement à droite sur {way_name:article}",
                "destination": "Tourner franchement à droite en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour",
                "name": "Faire demi-tour sur {way_name:article}",
                "destination": "Faire demi-tour en direction {destination:preposition}"
            }
        },
        "merge": {
            "default": {
                "default": "S’insérer {modifier}",
                "name": "S’insérer {modifier} sur {way_name:article}",
                "destination": "S’insérer {modifier} en direction {destination:preposition}"
            },
            "straight": {
                "default": "S’insérer",
                "name": "S’insérer sur {way_name:article}",
                "destination": "S’insérer en direction {destination:preposition}"
            },
            "slight left": {
                "default": "S’insérer légèrement à gauche",
                "name": "S’insérer légèrement à gauche sur {way_name:article}",
                "destination": "S’insérer légèrement à gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "S’insérer légèrement à droite",
                "name": "S’insérer légèrement à droite sur {way_name:article}",
                "destination": "S’insérer à droite en direction {destination:preposition}"
            },
            "sharp left": {
                "default": "S’insérer à gauche",
                "name": "S’insérer à gauche sur {way_name:article}",
                "destination": "S’insérer à gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "S’insérer à droite",
                "name": "S’insérer à droite sur {way_name:article}",
                "destination": "S’insérer à droite en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour",
                "name": "Faire demi-tour sur {way_name:article}",
                "destination": "Faire demi-tour en direction {destination:preposition}"
            }
        },
        "new name": {
            "default": {
                "default": "Continuer {modifier}",
                "name": "Continuer {modifier} sur {way_name:article}",
                "destination": "Continuer {modifier} en direction {destination:preposition}"
            },
            "straight": {
                "default": "Continuer tout droit",
                "name": "Continuer tout droit sur {way_name:article}",
                "destination": "Continuer tout droit en direction {destination:preposition}"
            },
            "sharp left": {
                "default": "Tourner franchement à gauche",
                "name": "Tourner franchement à gauche sur {way_name:article}",
                "destination": "Tourner franchement à gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "Tourner franchement à droite",
                "name": "Tourner franchement à droite sur {way_name:article}",
                "destination": "Tourner franchement à droite en direction {destination:preposition}"
            },
            "slight left": {
                "default": "Continuer légèrement à gauche",
                "name": "Continuer légèrement à gauche sur {way_name:article}",
                "destination": "Continuer légèrement à gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "Continuer légèrement à droite",
                "name": "Continuer légèrement à droite sur {way_name:article}",
                "destination": "Continuer légèrement à droite en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour",
                "name": "Faire demi-tour sur {way_name:article}",
                "destination": "Faire demi-tour en direction {destination:preposition}"
            }
        },
        "notification": {
            "default": {
                "default": "Continuer {modifier}",
                "name": "Continuer {modifier} sur {way_name:article}",
                "destination": "Continuer {modifier} en direction {destination:preposition}"
            },
            "uturn": {
                "default": "Faire demi-tour",
                "name": "Faire demi-tour sur {way_name:article}",
                "destination": "Faire demi-tour en direction {destination:preposition}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Prendre la sortie",
                "name": "Prendre la sortie sur {way_name:article}",
                "destination": "Prendre la sortie en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit}",
                "exit_destination": "Prendre la sortie {exit} en direction {destination:preposition}"
            },
            "left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la gauche",
                "exit_destination": "Prendre la sortie {exit} sur la gauche en direction {destination:preposition}"
            },
            "right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la droite",
                "exit_destination": "Prendre la sortie {exit} sur la droite en direction {destination:preposition}"
            },
            "sharp left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la gauche",
                "exit_destination": "Prendre la sortie {exit} sur la gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la droite",
                "exit_destination": "Prendre la sortie {exit} sur la droite en direction {destination:preposition}"
            },
            "slight left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la gauche",
                "exit_destination": "Prendre la sortie {exit} sur la gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}",
                "exit": "Prendre la sortie {exit} sur la droite",
                "exit_destination": "Prendre la sortie {exit} sur la droite en direction {destination:preposition}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Prendre la sortie",
                "name": "Prendre la sortie sur {way_name:article}",
                "destination": "Prendre la sortie en direction {destination:preposition}"
            },
            "left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}"
            },
            "right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}"
            },
            "sharp left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}"
            },
            "sharp right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}"
            },
            "slight left": {
                "default": "Prendre la sortie à gauche",
                "name": "Prendre la sortie à gauche sur {way_name:article}",
                "destination": "Prendre la sortie à gauche en direction {destination:preposition}"
            },
            "slight right": {
                "default": "Prendre la sortie à droite",
                "name": "Prendre la sortie à droite sur {way_name:article}",
                "destination": "Prendre la sortie à droite en direction {destination:preposition}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Prendre le rond-point",
                    "name": "Prendre le rond-point, puis sortir sur {way_name:article}",
                    "destination": "Prendre le rond-point, puis sortir en direction {destination:preposition}"
                },
                "name": {
                    "default": "Prendre {rotary_name:rotary}",
                    "name": "Prendre {rotary_name:rotary}, puis sortir par {way_name:article}",
                    "destination": "Prendre {rotary_name:rotary}, puis sortir en direction {destination:preposition}"
                },
                "exit": {
                    "default": "Prendre le rond-point, puis la {exit_number} sortie",
                    "name": "Prendre le rond-point, puis la {exit_number} sortie sur {way_name:article}",
                    "destination": "Prendre le rond-point, puis la {exit_number} sortie en direction {destination:preposition}"
                },
                "name_exit": {
                    "default": "Prendre {rotary_name:rotary}, puis la {exit_number} sortie",
                    "name": "Prendre {rotary_name:rotary}, puis la {exit_number} sortie sur {way_name:article}",
                    "destination": "Prendre {rotary_name:rotary}, puis la {exit_number} sortie en direction {destination:preposition}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Prendre le rond-point, puis la {exit_number} sortie",
                    "name": "Prendre le rond-point, puis la {exit_number} sortie sur {way_name:article}",
                    "destination": "Prendre le rond-point, puis la {exit_number} sortie en direction {destination:preposition}"
                },
                "default": {
                    "default": "Prendre le rond-point",
                    "name": "Prendre le rond-point, puis sortir sur {way_name:article}",
                    "destination": "Prendre le rond-point, puis sortir en direction {destination:preposition}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Tourner {modifier}",
                "name": "Tourner {modifier} sur {way_name:article}",
                "destination": "Tourner {modifier} en direction {destination:preposition}"
            },
            "left": {
                "default": "Tourner à gauche",
                "name": "Tourner à gauche sur {way_name:article}",
                "destination": "Tourner à gauche en direction {destination:preposition}"
            },
            "right": {
                "default": "Tourner à droite",
                "name": "Tourner à droite sur {way_name:article}",
                "destination": "Tourner à droite en direction {destination:preposition}"
            },
            "straight": {
                "default": "Continuer tout droit",
                "name": "Continuer tout droit sur {way_name:article}",
                "destination": "Continuer tout droit en direction {destination:preposition}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Sortir du rond-point",
                "name": "Sortir du rond-point sur {way_name:article}",
                "destination": "Sortir du rond-point en direction {destination:preposition}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Sortir du rond-point",
                "name": "Sortir du rond-point sur {way_name:article}",
                "destination": "Sortir du rond-point en direction {destination:preposition}"
            }
        },
        "turn": {
            "default": {
                "default": "Tourner {modifier}",
                "name": "Tourner {modifier} sur {way_name:article}",
                "destination": "Tourner {modifier} en direction {destination:preposition}"
            },
            "left": {
                "default": "Tourner à gauche",
                "name": "Tourner à gauche sur {way_name:article}",
                "destination": "Tourner à gauche en direction {destination:preposition}"
            },
            "right": {
                "default": "Tourner à droite",
                "name": "Tourner à droite sur {way_name:article}",
                "destination": "Tourner à droite en direction {destination:preposition}"
            },
            "straight": {
                "default": "Aller tout droit",
                "name": "Aller tout droit sur {way_name:article}",
                "destination": "Aller tout droit en direction {destination:preposition}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continuer tout droit"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],31:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "ראשונה",
                "2": "שניה",
                "3": "שלישית",
                "4": "רביעית",
                "5": "חמישית",
                "6": "שישית",
                "7": "שביעית",
                "8": "שמינית",
                "9": "תשיעית",
                "10": "עשירית"
            },
            "direction": {
                "north": "צפון",
                "northeast": "צפון מזרח",
                "east": "מזרח",
                "southeast": "דרום מזרח",
                "south": "דרום",
                "southwest": "דרום מערב",
                "west": "מערב",
                "northwest": "צפון מערב"
            },
            "modifier": {
                "left": "שמאלה",
                "right": "ימינה",
                "sharp left": "חדה שמאלה",
                "sharp right": "חדה ימינה",
                "slight left": "קלה שמאלה",
                "slight right": "קלה ימינה",
                "straight": "ישר",
                "uturn": "פניית פרסה"
            },
            "lanes": {
                "xo": "היצמד לימין",
                "ox": "היצמד לשמאל",
                "xox": "המשך בנתיב האמצעי",
                "oxo": "היצמד לימין או לשמאל"
            }
        },
        "modes": {
            "ferry": {
                "default": "עלה על המעבורת",
                "name": "עלה על המעבורת {way_name}",
                "destination": "עלה על המעבורת לכיוון {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, ואז, בעוד{distance}, {instruction_two}",
            "two linked": "{instruction_one}, ואז {instruction_two}",
            "one in distance": "בעוד {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "יציאה {exit}"
        },
        "arrive": {
            "default": {
                "default": "הגעת אל היעד ה{nth} שלך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name}"
            },
            "left": {
                "default": "הגעת אל היעד ה{nth} שלך משמאלך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך משמאלך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך משמאלך"
            },
            "right": {
                "default": "הגעת אל היעד ה{nth} שלך מימינך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך מימינך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך מימינך"
            },
            "sharp left": {
                "default": "הגעת אל היעד ה{nth} שלך משמאלך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך משמאלך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך משמאלך"
            },
            "sharp right": {
                "default": "הגעת אל היעד ה{nth} שלך מימינך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך מימינך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך מימינך"
            },
            "slight right": {
                "default": "הגעת אל היעד ה{nth} שלך מימינך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך מימינך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך מימינך"
            },
            "slight left": {
                "default": "הגעת אל היעד ה{nth} שלך משמאלך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך משמאלך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name} שלך משמאלך"
            },
            "straight": {
                "default": "הגעת אל היעד ה{nth} שלך, בהמשך",
                "upcoming": "אתה תגיע אל היעד ה{nth} שלך, בהמשך",
                "short": "הגעת",
                "short-upcoming": "תגיע",
                "named": "הגעת אל {waypoint_name}, בהמשך"
            }
        },
        "continue": {
            "default": {
                "default": "פנה {modifier}",
                "name": "פנה {modifier} כדי להישאר ב{way_name}",
                "destination": "פנה {modifier} לכיוון {destination}",
                "exit": "פנה {modifier} על {way_name}"
            },
            "straight": {
                "default": "המשך ישר",
                "name": "המשך ישר כדי להישאר על {way_name}",
                "destination": "המשך לכיוון {destination}",
                "distance": "המשך ישר לאורך {distance}",
                "namedistance": "המשך על {way_name} לאורך {distance}"
            },
            "sharp left": {
                "default": "פנה בחדות שמאלה",
                "name": "פנה בחדות שמאלה כדי להישאר על {way_name}",
                "destination": "פנה בחדות שמאלה לכיוון {destination}"
            },
            "sharp right": {
                "default": "פנה בחדות ימינה",
                "name": "פנה בחדות ימינה כדי להישאר על {way_name}",
                "destination": "פנה בחדות ימינה לכיוון {destination}"
            },
            "slight left": {
                "default": "פנה קלות שמאלה",
                "name": "פנה קלות שמאלה כדי להישאר על {way_name}",
                "destination": "פנה קלות שמאלה לכיוון {destination}"
            },
            "slight right": {
                "default": "פנה קלות ימינה",
                "name": "פנה קלות ימינה כדי להישאר על {way_name}",
                "destination": "פנה קלות ימינה לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה",
                "name": "פנה פניית פרסה והמשך על {way_name}",
                "destination": "פנה פניית פרסה לכיוון {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "התכוונן {direction}",
                "name": "התכוונן {direction} על {way_name}",
                "namedistance": "התכוונן {direction} על {way_name} לאורך {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "פנה {modifier}",
                "name": "פנה {modifier} על {way_name}",
                "destination": "פנה {modifier} לכיוון {destination}"
            },
            "straight": {
                "default": "המשך ישר",
                "name": "המשך ישר על {way_name}",
                "destination": "המשך ישר לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה בסוף הדרך",
                "name": "פנה פניית פרסה על {way_name} בסוף הדרך",
                "destination": "פנה פניית פרסה לכיוון {destination} בסוף הדרך"
            }
        },
        "fork": {
            "default": {
                "default": "היצמד {modifier} בהתפצלות",
                "name": "היצמד {modifier} על {way_name}",
                "destination": "היצמד {modifier} לכיוון {destination}"
            },
            "slight left": {
                "default": "היצמד לשמאל בהתפצלות",
                "name": "היצמד לשמאל על {way_name}",
                "destination": "היצמד לשמאל לכיוון {destination}"
            },
            "slight right": {
                "default": "היצמד ימינה בהתפצלות",
                "name": "היצמד לימין על {way_name}",
                "destination": "היצמד לימין לכיוון {destination}"
            },
            "sharp left": {
                "default": "פנה בחדות שמאלה בהתפצלות",
                "name": "פנה בחדות שמאלה על {way_name}",
                "destination": "פנה בחדות שמאלה לכיוון {destination}"
            },
            "sharp right": {
                "default": "פנה בחדות ימינה בהתפצלות",
                "name": "פנה בחדות ימינה על {way_name}",
                "destination": "פנה בחדות ימינה לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה",
                "name": "פנה פניית פרסה על {way_name}",
                "destination": "פנה פניית פרסה לכיוון {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "השתלב {modifier}",
                "name": "השתלב {modifier} על {way_name}",
                "destination": "השתלב {modifier} לכיוון {destination}"
            },
            "straight": {
                "default": "השתלב",
                "name": "השתלב על {way_name}",
                "destination": "השתלב לכיוון {destination}"
            },
            "slight left": {
                "default": "השתלב שמאלה",
                "name": "השתלב שמאלה על {way_name}",
                "destination": "השתלב שמאלה לכיוון {destination}"
            },
            "slight right": {
                "default": "השתלב ימינה",
                "name": "השתלב ימינה על {way_name}",
                "destination": "השתלב ימינה לכיוון {destination}"
            },
            "sharp left": {
                "default": "השתלב שמאלה",
                "name": "השתלב שמאלה על {way_name}",
                "destination": "השתלב שמאלה לכיוון {destination}"
            },
            "sharp right": {
                "default": "השתלב ימינה",
                "name": "השתלב ימינה על {way_name}",
                "destination": "השתלב ימינה לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה",
                "name": "פנה פניית פרסה על {way_name}",
                "destination": "פנה פניית פרסה לכיוון {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "המשך {modifier}",
                "name": "המשך {modifier} על {way_name}",
                "destination": "המשך {modifier} לכיוון {destination}"
            },
            "straight": {
                "default": "המשך ישר",
                "name": "המשך על {way_name}",
                "destination": "המשך לכיוון {destination}"
            },
            "sharp left": {
                "default": "פנה בחדות שמאלה",
                "name": "פנה בחדות שמאלה על {way_name}",
                "destination": "פנה בחדות שמאלה לכיוון {destination}"
            },
            "sharp right": {
                "default": "פנה בחדות ימינה",
                "name": "פנה בחדות ימינה על {way_name}",
                "destination": "פנה בחדות ימינה לכיוון {destination}"
            },
            "slight left": {
                "default": "המשך בנטייה קלה שמאלה",
                "name": "המשך בנטייה קלה שמאלה על {way_name}",
                "destination": "המשך בנטייה קלה שמאלה לכיוון {destination}"
            },
            "slight right": {
                "default": "המשך בנטייה קלה ימינה",
                "name": "המשך בנטייה קלה ימינה על {way_name}",
                "destination": "המשך בנטייה קלה ימינה לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה",
                "name": "פנה פניית פרסה על {way_name}",
                "destination": "פנה פניית פרסה לכיוון {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "המשך {modifier}",
                "name": "המשך {modifier} על {way_name}",
                "destination": "המשך {modifier} לכיוון {destination}"
            },
            "uturn": {
                "default": "פנה פניית פרסה",
                "name": "פנה פניית פרסה על {way_name}",
                "destination": "פנה פניית פרסה לכיוון {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "צא ביציאה",
                "name": "צא ביציאה על {way_name}",
                "destination": "צא ביציאה לכיוון {destination}",
                "exit": "צא ביציאה {exit}",
                "exit_destination": "צא ביציאה {exit} לכיוון {destination}"
            },
            "left": {
                "default": "צא ביציאה שמשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}",
                "exit": "צא ביציאה {exit} משמאלך",
                "exit_destination": "צא ביציאה {exit} משמאלך לכיוון {destination}"
            },
            "right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}",
                "exit": "צא ביציאה {exit} מימינך",
                "exit_destination": "צא ביציאה {exit} מימינך לכיוון {destination}"
            },
            "sharp left": {
                "default": "צא ביציאה שבשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}",
                "exit": "צא ביציאה {exit} משמאלך",
                "exit_destination": "צא ביציאה {exit} משמאלך לכיוון {destination}"
            },
            "sharp right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}",
                "exit": "צא ביציאה {exit} מימינך",
                "exit_destination": "צא ביציאה {exit} מימינך לכיוון {destination}"
            },
            "slight left": {
                "default": "צא ביציאה שבשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}",
                "exit": "צא ביציאה {exit} משמאלך",
                "exit_destination": "צא ביציאה {exit} משמאלך לכיוון {destination}"
            },
            "slight right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}",
                "exit": "צא ביציאה {exit} מימינך",
                "exit_destination": "צא ביציאה {exit} מימינך לכיוון {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "צא ביציאה",
                "name": "צא ביציאה על {way_name}",
                "destination": "צא ביציאה לכיוון {destination}"
            },
            "left": {
                "default": "צא ביציאה שבשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}"
            },
            "right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}"
            },
            "sharp left": {
                "default": "צא ביציאה שבשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}"
            },
            "sharp right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}"
            },
            "slight left": {
                "default": "צא ביציאה שבשמאלך",
                "name": "צא ביציאה שמשמאלך על {way_name}",
                "destination": "צא ביציאה שמשמאלך לכיוון {destination}"
            },
            "slight right": {
                "default": "צא ביציאה שמימינך",
                "name": "צא ביציאה שמימינך על {way_name}",
                "destination": "צא ביציאה שמימינך לכיוון {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "השתלב במעגל התנועה",
                    "name": "השתלב במעגל התנועה וצא על {way_name}",
                    "destination": "השתלב במעגל התנועה וצא לכיוון {destination}"
                },
                "name": {
                    "default": "היכנס ל{rotary_name}",
                    "name": "היכנס ל{rotary_name} וצא על {way_name}",
                    "destination": "היכנס ל{rotary_name} וצא לכיוון {destination}"
                },
                "exit": {
                    "default": "השתלב במעגל התנועה וצא ביציאה {exit_number}",
                    "name": "השתלב במעגל התנועה וצא ביציאה {exit_number} ל{way_name}",
                    "destination": "השתלב במעגל התנועה וצא ביציאה {exit_number} לכיוון {destination}"
                },
                "name_exit": {
                    "default": "היכנס ל{rotary_name} וצא ביציאה ה{exit_number}",
                    "name": "היכנס ל{rotary_name} וצא ביציאה ה{exit_number} ל{way_name}",
                    "destination": "היכנס ל{rotary_name} וצא ביציאה ה{exit_number} לכיוון {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "השתלב במעגל התנועה וצא ביציאה {exit_number}",
                    "name": "השתלב במעגל התנועה וצא ביציאה {exit_number} ל{way_name}",
                    "destination": "השתלב במעגל התנועה וצא ביציאה {exit_number} לכיוון {destination}"
                },
                "default": {
                    "default": "השתלב במעגל התנועה",
                    "name": "השתלב במעגל התנועה וצא על {way_name}",
                    "destination": "השתלב במעגל התנועה וצא לכיוון {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "פנה {modifier}",
                "name": "פנה {modifier} על {way_name}",
                "destination": "פנה {modifier} לכיוון {destination}"
            },
            "left": {
                "default": "פנה שמאלה",
                "name": "פנה שמאלה ל{way_name}",
                "destination": "פנה שמאלה לכיוון {destination}"
            },
            "right": {
                "default": "פנה ימינה",
                "name": "פנה ימינה ל{way_name}",
                "destination": "פנה ימינה לכיוון {destination}"
            },
            "straight": {
                "default": "המשך ישר",
                "name": "המשך ישר על {way_name}",
                "destination": "המשך ישר לכיוון {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "צא ממעגל התנועה",
                "name": "צא ממעגל התנועה ל{way_name}",
                "destination": "צא ממעגל התנועה לכיוון {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "צא ממעגל התנועה",
                "name": "צא ממעגל התנועה ל{way_name}",
                "destination": "צא ממעגל התנועה לכיוון {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "פנה {modifier}",
                "name": "פנה {modifier} על {way_name}",
                "destination": "פנה {modifier} לכיוון {destination}"
            },
            "left": {
                "default": "פנה שמאלה",
                "name": "פנה שמאלה ל{way_name}",
                "destination": "פנה שמאלה לכיוון {destination}"
            },
            "right": {
                "default": "פנה ימינה",
                "name": "פנה ימינה ל{way_name}",
                "destination": "פנה ימינה לכיוון {destination}"
            },
            "straight": {
                "default": "המשך ישר",
                "name": "המשך ישר ל{way_name}",
                "destination": "המשך ישר לכיוון {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "המשך ישר"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],32:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1",
                "2": "2",
                "3": "3",
                "4": "4",
                "5": "5",
                "6": "6",
                "7": "7",
                "8": "8",
                "9": "9",
                "10": "10"
            },
            "direction": {
                "north": "utara",
                "northeast": "timur laut",
                "east": "timur",
                "southeast": "tenggara",
                "south": "selatan",
                "southwest": "barat daya",
                "west": "barat",
                "northwest": "barat laut"
            },
            "modifier": {
                "left": "kiri",
                "right": "kanan",
                "sharp left": "tajam kiri",
                "sharp right": "tajam kanan",
                "slight left": "agak ke kiri",
                "slight right": "agak ke kanan",
                "straight": "lurus",
                "uturn": "putar balik"
            },
            "lanes": {
                "xo": "Tetap di kanan",
                "ox": "Tetap di kiri",
                "xox": "Tetap di tengah",
                "oxo": "Tetap di kiri atau kanan"
            }
        },
        "modes": {
            "ferry": {
                "default": "Naik ferry",
                "name": "Naik ferry di {way_name}",
                "destination": "Naik ferry menuju {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, then, in {distance}, {instruction_two}",
            "two linked": "{instruction_one}, then {instruction_two}",
            "one in distance": "In {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "Anda telah tiba di tujuan ke-{nth}",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}"
            },
            "left": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kiri"
            },
            "right": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kanan"
            },
            "sharp left": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kiri"
            },
            "sharp right": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kanan"
            },
            "slight right": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kanan",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kanan"
            },
            "slight left": {
                "default": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, di sebelah kiri",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, di sebelah kiri"
            },
            "straight": {
                "default": "Anda telah tiba di tujuan ke-{nth}, lurus saja",
                "upcoming": "Anda telah tiba di tujuan ke-{nth}, lurus saja",
                "short": "Anda telah tiba di tujuan ke-{nth}",
                "short-upcoming": "Anda telah tiba di tujuan ke-{nth}",
                "named": "Anda telah tiba di {waypoint_name}, lurus saja"
            }
        },
        "continue": {
            "default": {
                "default": "Belok {modifier}",
                "name": "Terus {modifier} ke {way_name}",
                "destination": "Belok {modifier} menuju {destination}",
                "exit": "Belok {modifier} ke {way_name}"
            },
            "straight": {
                "default": "Lurus terus",
                "name": "Terus ke {way_name}",
                "destination": "Terus menuju {destination}",
                "distance": "Continue straight for {distance}",
                "namedistance": "Continue on {way_name} for {distance}"
            },
            "sharp left": {
                "default": "Belok kiri tajam",
                "name": "Make a sharp left to stay on {way_name}",
                "destination": "Belok kiri tajam menuju {destination}"
            },
            "sharp right": {
                "default": "Belok kanan tajam",
                "name": "Make a sharp right to stay on {way_name}",
                "destination": "Belok kanan tajam menuju {destination}"
            },
            "slight left": {
                "default": "Tetap agak di kiri",
                "name": "Tetap agak di kiri ke {way_name}",
                "destination": "Tetap agak di kiri menuju {destination}"
            },
            "slight right": {
                "default": "Tetap agak di kanan",
                "name": "Tetap agak di kanan ke {way_name}",
                "destination": "Tetap agak di kanan menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik",
                "name": "Putar balik ke arah {way_name}",
                "destination": "Putar balik menuju {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Arah {direction}",
                "name": "Arah {direction} di {way_name}",
                "namedistance": "Head {direction} on {way_name} for {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Belok {modifier}",
                "name": "Belok {modifier} ke {way_name}",
                "destination": "Belok {modifier} menuju {destination}"
            },
            "straight": {
                "default": "Lurus terus",
                "name": "Tetap lurus ke {way_name} ",
                "destination": "Tetap lurus menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik di akhir jalan",
                "name": "Putar balik di {way_name} di akhir jalan",
                "destination": "Putar balik menuju {destination} di akhir jalan"
            }
        },
        "fork": {
            "default": {
                "default": "Tetap {modifier} di pertigaan",
                "name": "Tetap {modifier} di pertigaan ke {way_name}",
                "destination": "Tetap {modifier} di pertigaan menuju {destination}"
            },
            "slight left": {
                "default": "Tetap di kiri pada pertigaan",
                "name": "Tetap di kiri pada pertigaan ke arah {way_name}",
                "destination": "Tetap di kiri pada pertigaan menuju {destination}"
            },
            "slight right": {
                "default": "Tetap di kanan pada pertigaan",
                "name": "Tetap di kanan pada pertigaan ke arah {way_name}",
                "destination": "Tetap di kanan pada pertigaan menuju {destination}"
            },
            "sharp left": {
                "default": "Belok kiri pada pertigaan",
                "name": "Belok kiri tajam ke arah {way_name}",
                "destination": "Belok kiri tajam menuju {destination}"
            },
            "sharp right": {
                "default": "Belok kanan pada pertigaan",
                "name": "Belok kanan tajam ke arah {way_name}",
                "destination": "Belok kanan tajam menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik",
                "name": "Putar balik ke arah {way_name}",
                "destination": "Putar balik menuju {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Bergabung {modifier}",
                "name": "Bergabung {modifier} ke arah {way_name}",
                "destination": "Bergabung {modifier} menuju {destination}"
            },
            "straight": {
                "default": "Bergabung lurus",
                "name": "Bergabung lurus ke arah {way_name}",
                "destination": "Bergabung lurus menuju {destination}"
            },
            "slight left": {
                "default": "Bergabung di kiri",
                "name": "Bergabung di kiri ke arah {way_name}",
                "destination": "Bergabung di kiri menuju {destination}"
            },
            "slight right": {
                "default": "Bergabung di kanan",
                "name": "Bergabung di kanan ke arah {way_name}",
                "destination": "Bergabung di kanan menuju {destination}"
            },
            "sharp left": {
                "default": "Bergabung di kiri",
                "name": "Bergabung di kiri ke arah {way_name}",
                "destination": "Bergabung di kiri menuju {destination}"
            },
            "sharp right": {
                "default": "Bergabung di kanan",
                "name": "Bergabung di kanan ke arah {way_name}",
                "destination": "Bergabung di kanan menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik",
                "name": "Putar balik ke arah {way_name}",
                "destination": "Putar balik menuju {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Lanjutkan {modifier}",
                "name": "Lanjutkan {modifier} menuju {way_name}",
                "destination": "Lanjutkan {modifier} menuju {destination}"
            },
            "straight": {
                "default": "Lurus terus",
                "name": "Terus ke {way_name}",
                "destination": "Terus menuju {destination}"
            },
            "sharp left": {
                "default": "Belok kiri tajam",
                "name": "Belok kiri tajam ke arah {way_name}",
                "destination": "Belok kiri tajam menuju {destination}"
            },
            "sharp right": {
                "default": "Belok kanan tajam",
                "name": "Belok kanan tajam ke arah {way_name}",
                "destination": "Belok kanan tajam menuju {destination}"
            },
            "slight left": {
                "default": "Lanjut dengan agak ke kiri",
                "name": "Lanjut dengan agak di kiri ke {way_name}",
                "destination": "Tetap agak di kiri menuju {destination}"
            },
            "slight right": {
                "default": "Tetap agak di kanan",
                "name": "Tetap agak di kanan ke {way_name}",
                "destination": "Tetap agak di kanan menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik",
                "name": "Putar balik ke arah {way_name}",
                "destination": "Putar balik menuju {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Lanjutkan {modifier}",
                "name": "Lanjutkan {modifier} menuju {way_name}",
                "destination": "Lanjutkan {modifier} menuju {destination}"
            },
            "uturn": {
                "default": "Putar balik",
                "name": "Putar balik ke arah {way_name}",
                "destination": "Putar balik menuju {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Ambil jalan melandai",
                "name": "Ambil jalan melandai ke {way_name}",
                "destination": "Ambil jalan melandai menuju {destination}",
                "exit": "Take exit {exit}",
                "exit_destination": "Take exit {exit} towards {destination}"
            },
            "left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan menuju {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            },
            "sharp left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "sharp right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan menuju {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            },
            "slight left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}",
                "exit": "Take exit {exit} on the left",
                "exit_destination": "Take exit {exit} on the left towards {destination}"
            },
            "slight right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan  menuju {destination}",
                "exit": "Take exit {exit} on the right",
                "exit_destination": "Take exit {exit} on the right towards {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Ambil jalan melandai",
                "name": "Ambil jalan melandai ke {way_name}",
                "destination": "Ambil jalan melandai menuju {destination}"
            },
            "left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}"
            },
            "right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan  menuju {destination}"
            },
            "sharp left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}"
            },
            "sharp right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan  menuju {destination}"
            },
            "slight left": {
                "default": "Ambil jalan yang melandai di sebelah kiri",
                "name": "Ambil jalan melandai di sebelah kiri ke arah {way_name}",
                "destination": "Ambil jalan melandai di sebelah kiri menuju {destination}"
            },
            "slight right": {
                "default": "Ambil jalan melandai di sebelah kanan",
                "name": "Ambil jalan melandai di sebelah kanan ke {way_name}",
                "destination": "Ambil jalan melandai di sebelah kanan  menuju {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Masuk bundaran",
                    "name": "Masuk bundaran dan keluar arah {way_name}",
                    "destination": "Masuk bundaran dan keluar menuju {destination}"
                },
                "name": {
                    "default": "Masuk {rotary_name}",
                    "name": "Masuk {rotary_name} dan keluar arah {way_name}",
                    "destination": "Masuk {rotary_name} dan keluar menuju {destination}"
                },
                "exit": {
                    "default": "Masuk bundaran dan ambil jalan keluar {exit_number}",
                    "name": "Masuk bundaran dan ambil jalan keluar {exit_number} arah {way_name}",
                    "destination": "Masuk bundaran dan ambil jalan keluar {exit_number} menuju {destination}"
                },
                "name_exit": {
                    "default": "Masuk {rotary_name} dan ambil jalan keluar {exit_number}",
                    "name": "Masuk {rotary_name} dan ambil jalan keluar {exit_number} arah {way_name}",
                    "destination": "Masuk {rotary_name} dan ambil jalan keluar {exit_number} menuju {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Masuk bundaran dan ambil jalan keluar {exit_number}",
                    "name": "Masuk bundaran dan ambil jalan keluar {exit_number} arah {way_name}",
                    "destination": "Masuk bundaran dan ambil jalan keluar {exit_number} menuju {destination}"
                },
                "default": {
                    "default": "Masuk bundaran",
                    "name": "Masuk bundaran dan keluar arah {way_name}",
                    "destination": "Masuk bundaran dan keluar menuju {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Lakukan {modifier}",
                "name": "Lakukan {modifier} ke arah {way_name}",
                "destination": "Lakukan {modifier} menuju {destination}"
            },
            "left": {
                "default": "Belok kiri",
                "name": "Belok kiri ke {way_name}",
                "destination": "Belok kiri menuju {destination}"
            },
            "right": {
                "default": "Belok kanan",
                "name": "Belok kanan ke {way_name}",
                "destination": "Belok kanan menuju {destination}"
            },
            "straight": {
                "default": "Lurus terus",
                "name": "Tetap lurus ke {way_name} ",
                "destination": "Tetap lurus menuju {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Lakukan {modifier}",
                "name": "Lakukan {modifier} ke arah {way_name}",
                "destination": "Lakukan {modifier} menuju {destination}"
            },
            "left": {
                "default": "Belok kiri",
                "name": "Belok kiri ke {way_name}",
                "destination": "Belok kiri menuju {destination}"
            },
            "right": {
                "default": "Belok kanan",
                "name": "Belok kanan ke {way_name}",
                "destination": "Belok kanan menuju {destination}"
            },
            "straight": {
                "default": "Lurus terus",
                "name": "Tetap lurus ke {way_name} ",
                "destination": "Tetap lurus menuju {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Lakukan {modifier}",
                "name": "Lakukan {modifier} ke arah {way_name}",
                "destination": "Lakukan {modifier} menuju {destination}"
            },
            "left": {
                "default": "Belok kiri",
                "name": "Belok kiri ke {way_name}",
                "destination": "Belok kiri menuju {destination}"
            },
            "right": {
                "default": "Belok kanan",
                "name": "Belok kanan ke {way_name}",
                "destination": "Belok kanan menuju {destination}"
            },
            "straight": {
                "default": "Lurus",
                "name": "Lurus arah {way_name}",
                "destination": "Lurus menuju {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Lakukan {modifier}",
                "name": "Lakukan {modifier} ke arah {way_name}",
                "destination": "Lakukan {modifier} menuju {destination}"
            },
            "left": {
                "default": "Belok kiri",
                "name": "Belok kiri ke {way_name}",
                "destination": "Belok kiri menuju {destination}"
            },
            "right": {
                "default": "Belok kanan",
                "name": "Belok kanan ke {way_name}",
                "destination": "Belok kanan menuju {destination}"
            },
            "straight": {
                "default": "Lurus",
                "name": "Lurus arah {way_name}",
                "destination": "Lurus menuju {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Lurus terus"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],33:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1ª",
                "2": "2ª",
                "3": "3ª",
                "4": "4ª",
                "5": "5ª",
                "6": "6ª",
                "7": "7ª",
                "8": "8ª",
                "9": "9ª",
                "10": "10ª"
            },
            "direction": {
                "north": "nord",
                "northeast": "nord-est",
                "east": "est",
                "southeast": "sud-est",
                "south": "sud",
                "southwest": "sud-ovest",
                "west": "ovest",
                "northwest": "nord-ovest"
            },
            "modifier": {
                "left": "sinistra",
                "right": "destra",
                "sharp left": "sinistra",
                "sharp right": "destra",
                "slight left": "sinistra leggermente",
                "slight right": "destra leggermente",
                "straight": "dritto",
                "uturn": "inversione a U"
            },
            "lanes": {
                "xo": "Mantieni la destra",
                "ox": "Mantieni la sinistra",
                "xox": "Rimani in mezzo",
                "oxo": "Mantieni la destra o la sinistra"
            }
        },
        "modes": {
            "ferry": {
                "default": "Prendi il traghetto",
                "name": "Prendi il traghetto {way_name}",
                "destination": "Prendi il traghetto verso {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, poi tra {distance},{instruction_two}",
            "two linked": "{instruction_one}, poi {instruction_two}",
            "one in distance": "tra {distance} {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "Sei arrivato alla tua {nth} destinazione",
                "upcoming": "Sei arrivato alla tua {nth} destinazione",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "Sei arrivato a {waypoint_name}"
            },
            "left": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla sinistra"
            },
            "right": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla destra"
            },
            "sharp left": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla sinistra"
            },
            "sharp right": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla destra"
            },
            "slight right": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla destra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla destra"
            },
            "slight left": {
                "default": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "upcoming": "sei arrivato alla tua {nth} destinazione, sulla sinistra",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, sulla sinistra"
            },
            "straight": {
                "default": "sei arrivato alla tua {nth} destinazione, si trova davanti a te",
                "upcoming": "sei arrivato alla tua {nth} destinazione, si trova davanti a te",
                "short": "Sei arrivato alla tua {nth} destinazione",
                "short-upcoming": "Sei arrivato alla tua {nth} destinazione",
                "named": "sei arrivato a {waypoint_name}, si trova davanti a te"
            }
        },
        "continue": {
            "default": {
                "default": "Gira a {modifier}",
                "name": "Gira a {modifier} per stare su {way_name}",
                "destination": "Gira a {modifier} verso {destination}",
                "exit": "Gira a {modifier} in {way_name}"
            },
            "straight": {
                "default": "Continua dritto",
                "name": "Continua dritto per stare su {way_name}",
                "destination": "Continua verso {destination}",
                "distance": "Continua dritto per {distance}",
                "namedistance": "Continua su {way_name} per {distance}"
            },
            "sharp left": {
                "default": "Svolta a sinistra",
                "name": "Fai una stretta curva a sinistra per stare su {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Svolta a destra",
                "name": "Fau una stretta curva a destra per stare su {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "slight left": {
                "default": "Fai una leggera curva a sinistra",
                "name": "Fai una leggera curva a sinistra per stare su {way_name}",
                "destination": "Fai una leggera curva a sinistra verso {destination}"
            },
            "slight right": {
                "default": "Fai una leggera curva a destra",
                "name": "Fai una leggera curva a destra per stare su {way_name}",
                "destination": "Fai una leggera curva a destra verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U",
                "name": "Fai un'inversione ad U poi continua su {way_name}",
                "destination": "Fai un'inversione a U verso {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Continua verso {direction}",
                "name": "Continua verso {direction} in {way_name}",
                "namedistance": "Head {direction} on {way_name} for {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Gira a {modifier}",
                "name": "Gira a {modifier} in {way_name}",
                "destination": "Gira a {modifier} verso {destination}"
            },
            "straight": {
                "default": "Continua dritto",
                "name": "Continua dritto in {way_name}",
                "destination": "Continua dritto verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U alla fine della strada",
                "name": "Fai un'inversione a U in {way_name} alla fine della strada",
                "destination": "Fai un'inversione a U verso {destination} alla fine della strada"
            }
        },
        "fork": {
            "default": {
                "default": "Mantieni la {modifier} al bivio",
                "name": "Mantieni la {modifier} al bivio in {way_name}",
                "destination": "Mantieni la {modifier} al bivio verso {destination}"
            },
            "slight left": {
                "default": "Mantieni la sinistra al bivio",
                "name": "Mantieni la sinistra al bivio in {way_name}",
                "destination": "Mantieni la sinistra al bivio verso {destination}"
            },
            "slight right": {
                "default": "Mantieni la destra al bivio",
                "name": "Mantieni la destra al bivio in {way_name}",
                "destination": "Mantieni la destra al bivio verso {destination}"
            },
            "sharp left": {
                "default": "Svolta a sinistra al bivio",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Svolta a destra al bivio",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U",
                "name": "Fai un'inversione a U in {way_name}",
                "destination": "Fai un'inversione a U verso {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Immettiti a {modifier}",
                "name": "Immettiti {modifier} in {way_name}",
                "destination": "Immettiti {modifier} verso {destination}"
            },
            "straight": {
                "default": "Immettiti a dritto",
                "name": "Immettiti dritto in {way_name}",
                "destination": "Immettiti dritto verso {destination}"
            },
            "slight left": {
                "default": "Immettiti a sinistra",
                "name": "Immettiti a sinistra in {way_name}",
                "destination": "Immettiti a sinistra verso {destination}"
            },
            "slight right": {
                "default": "Immettiti a destra",
                "name": "Immettiti a destra in {way_name}",
                "destination": "Immettiti a destra verso {destination}"
            },
            "sharp left": {
                "default": "Immettiti a sinistra",
                "name": "Immettiti a sinistra in {way_name}",
                "destination": "Immettiti a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Immettiti a destra",
                "name": "Immettiti a destra in {way_name}",
                "destination": "Immettiti a destra verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U",
                "name": "Fai un'inversione a U in {way_name}",
                "destination": "Fai un'inversione a U verso {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continua a {modifier}",
                "name": "Continua a {modifier} in {way_name}",
                "destination": "Continua a {modifier} verso {destination}"
            },
            "straight": {
                "default": "Continua dritto",
                "name": "Continua in {way_name}",
                "destination": "Continua verso {destination}"
            },
            "sharp left": {
                "default": "Svolta a sinistra",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Svolta a destra",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "slight left": {
                "default": "Continua leggermente a sinistra",
                "name": "Continua leggermente a sinistra in {way_name}",
                "destination": "Continua leggermente a sinistra verso {destination}"
            },
            "slight right": {
                "default": "Continua leggermente a destra",
                "name": "Continua leggermente a destra in {way_name} ",
                "destination": "Continua leggermente a destra verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U",
                "name": "Fai un'inversione a U in {way_name}",
                "destination": "Fai un'inversione a U verso {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continua a {modifier}",
                "name": "Continua a {modifier} in {way_name}",
                "destination": "Continua a {modifier} verso {destination}"
            },
            "uturn": {
                "default": "Fai un'inversione a U",
                "name": "Fai un'inversione a U in {way_name}",
                "destination": "Fai un'inversione a U verso {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Prendi la rampa",
                "name": "Prendi la rampa in {way_name}",
                "destination": "Prendi la rampa verso {destination}",
                "exit": "Prendi l'uscita {exit}",
                "exit_destination": "Prendi l'uscita  {exit} verso {destination}"
            },
            "left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}",
                "exit": "Prendi l'uscita {exit} a sinistra",
                "exit_destination": "Prendi la {exit}  uscita a sinistra verso {destination}"
            },
            "right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}",
                "exit": "Prendi la {exit} uscita a destra",
                "exit_destination": "Prendi la {exit} uscita a destra verso {destination}"
            },
            "sharp left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}",
                "exit": "Prendi l'uscita {exit} a sinistra",
                "exit_destination": "Prendi la {exit}  uscita a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}",
                "exit": "Prendi la {exit} uscita a destra",
                "exit_destination": "Prendi la {exit} uscita a destra verso {destination}"
            },
            "slight left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}",
                "exit": "Prendi l'uscita {exit} a sinistra",
                "exit_destination": "Prendi la {exit}  uscita a sinistra verso {destination}"
            },
            "slight right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}",
                "exit": "Prendi la {exit} uscita a destra",
                "exit_destination": "Prendi la {exit} uscita a destra verso {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Prendi la rampa",
                "name": "Prendi la rampa in {way_name}",
                "destination": "Prendi la rampa verso {destination}"
            },
            "left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}"
            },
            "right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}"
            },
            "sharp left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}"
            },
            "sharp right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}"
            },
            "slight left": {
                "default": "Prendi la rampa a sinistra",
                "name": "Prendi la rampa a sinistra in {way_name}",
                "destination": "Prendi la rampa a sinistra verso {destination}"
            },
            "slight right": {
                "default": "Prendi la rampa a destra",
                "name": "Prendi la rampa a destra in {way_name}",
                "destination": "Prendi la rampa a destra verso {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Immettiti nella rotonda",
                    "name": "Immettiti nella ritonda ed esci in {way_name}",
                    "destination": "Immettiti nella ritonda ed esci verso {destination}"
                },
                "name": {
                    "default": "Immettiti in {rotary_name}",
                    "name": "Immettiti in {rotary_name} ed esci su {way_name}",
                    "destination": "Immettiti in {rotary_name} ed esci verso {destination}"
                },
                "exit": {
                    "default": "Immettiti nella rotonda e prendi la {exit_number} uscita",
                    "name": "Immettiti nella rotonda e prendi la {exit_number} uscita in {way_name}",
                    "destination": "Immettiti nella rotonda e prendi la {exit_number} uscita verso   {destination}"
                },
                "name_exit": {
                    "default": "Immettiti in {rotary_name} e prendi la {exit_number} uscita",
                    "name": "Immettiti in {rotary_name} e prendi la {exit_number} uscita in {way_name}",
                    "destination": "Immettiti in {rotary_name} e prendi la {exit_number}  uscita verso {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Immettiti nella rotonda e prendi la {exit_number} uscita",
                    "name": "Immettiti nella rotonda e prendi la {exit_number} uscita in {way_name}",
                    "destination": "Immettiti nella rotonda e prendi la {exit_number} uscita verso {destination}"
                },
                "default": {
                    "default": "Entra nella rotonda",
                    "name": "Entra nella rotonda e prendi l'uscita in {way_name}",
                    "destination": "Entra nella rotonda e prendi l'uscita verso {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Fai una {modifier}",
                "name": "Fai una {modifier} in {way_name}",
                "destination": "Fai una {modifier} verso {destination}"
            },
            "left": {
                "default": "Svolta a sinistra",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "right": {
                "default": "Gira a destra",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "straight": {
                "default": "Continua dritto",
                "name": "Continua dritto in {way_name}",
                "destination": "Continua dritto verso {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Fai una {modifier}",
                "name": "Fai una {modifier} in {way_name}",
                "destination": "Fai una {modifier} verso {destination}"
            },
            "left": {
                "default": "Svolta a sinistra",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "right": {
                "default": "Gira a destra",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "straight": {
                "default": "Continua dritto",
                "name": "Continua dritto in {way_name}",
                "destination": "Continua dritto verso {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Fai una {modifier}",
                "name": "Fai una {modifier} in {way_name}",
                "destination": "Fai una {modifier} verso {destination}"
            },
            "left": {
                "default": "Svolta a sinistra",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "right": {
                "default": "Gira a destra",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "straight": {
                "default": "Prosegui dritto",
                "name": "Continua su {way_name}",
                "destination": "Continua verso {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Fai una {modifier}",
                "name": "Fai una {modifier} in {way_name}",
                "destination": "Fai una {modifier} verso {destination}"
            },
            "left": {
                "default": "Svolta a sinistra",
                "name": "Svolta a sinistra in {way_name}",
                "destination": "Svolta a sinistra verso {destination}"
            },
            "right": {
                "default": "Gira a destra",
                "name": "Svolta a destra in {way_name}",
                "destination": "Svolta a destra verso {destination}"
            },
            "straight": {
                "default": "Prosegui dritto",
                "name": "Continua su {way_name}",
                "destination": "Continua verso {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continua dritto"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],34:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": false
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "첫번쩨",
                "2": "두번째",
                "3": "세번째",
                "4": "네번쩨",
                "5": "다섯번째",
                "6": "여섯번째",
                "7": "일곱번째",
                "8": "여덟번째",
                "9": "아홉번째",
                "10": "열번째"
            },
            "direction": {
                "north": "북쪽",
                "northeast": "북동쪽",
                "east": "동쪽",
                "southeast": "남동쪽",
                "south": "남쪽",
                "southwest": "남서쪽",
                "west": "서쪽",
                "northwest": "북서쪽"
            },
            "modifier": {
                "left": "좌회전",
                "right": "우회전",
                "sharp left": "바로좌회전",
                "sharp right": "바로우회전",
                "slight left": "조금왼쪽",
                "slight right": "조금오른쪽",
                "straight": "직진",
                "uturn": "유턴"
            },
            "lanes": {
                "xo": "우측차선 유지",
                "ox": "좌측차선 유지",
                "xox": "중앙유지",
                "oxo": "계속 좌측 또는 우측 차선"
            }
        },
        "modes": {
            "ferry": {
                "default": "페리를 타시오",
                "name": "페리를 타시오 {way_name}",
                "destination": "페리를 타고 {destination}까지 가세요."
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, 그리고, {distance} 안에, {instruction_two}",
            "two linked": "{instruction_one}, 그리고 {instruction_two}",
            "one in distance": "{distance} 내에, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "{exit}번으로 나가세요."
        },
        "arrive": {
            "default": {
                "default": " {nth}목적지에 도착하였습니다.",
                "upcoming": "{nth}목적지에 곧 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "도착할 예정입니다.",
                "named": "경유지 {waypoint_name}에 도착하였습니다."
            },
            "left": {
                "default": "좌측에 {nth} 목적지가 있습니다.",
                "upcoming": "좌측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "좌측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "right": {
                "default": "우측에 {nth} 목적지가 있습니다.",
                "upcoming": "우측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "우측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "sharp left": {
                "default": "좌측에 {nth} 목적지가 있습니다.",
                "upcoming": "좌측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "좌측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "sharp right": {
                "default": "우측에 {nth} 목적지가 있습니다.",
                "upcoming": "우측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "우측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "slight right": {
                "default": "우측에 {nth} 목적지가 있습니다.",
                "upcoming": "우측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "우측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "slight left": {
                "default": "좌측에 {nth} 목적지가 있습니다.",
                "upcoming": "좌측에 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "좌측에 경유지 {waypoint_name}에 도착하였습니다."
            },
            "straight": {
                "default": "바로 앞에 {nth} 목적지가 있습니다.",
                "upcoming": "직진하시면 {nth} 목적지에 도착할 예정입니다.",
                "short": "도착하였습니다",
                "short-upcoming": "목적지에 곧 도착할 예정입니다.",
                "named": "정면에 경유지 {waypoint_name}에 도착하였습니다."
            }
        },
        "continue": {
            "default": {
                "default": "{modifier} 회전",
                "name": "{modifier} 회전하고 {way_name}로 직진해 주세요.",
                "destination": "{modifier} 회전하고 {destination}까지 가세요.",
                "exit": "{way_name} 쪽으로 {modifier} 회전 하세요."
            },
            "straight": {
                "default": "계속 직진해 주세요.",
                "name": "{way_name} 로 계속 직진해 주세요.",
                "destination": "{destination}까지 직진해 주세요.",
                "distance": "{distance}까지 직진해 주세요.",
                "namedistance": "{distance}까지 {way_name}로 가주세요."
            },
            "sharp left": {
                "default": "급좌회전 하세요.",
                "name": "급좌회전 하신 후 {way_name}로 가세요.",
                "destination": "급좌회전 하신 후 {destination}로 가세요."
            },
            "sharp right": {
                "default": "급우회전 하세요.",
                "name": "급우회전 하고 {way_name}로 가세요.",
                "destination": "급우회전 하신 후 {destination}로 가세요."
            },
            "slight left": {
                "default": "약간 좌회전하세요.",
                "name": "약간 좌회전 하고 {way_name}로 가세요.",
                "destination": "약간 좌회전 하신 후 {destination}로 가세요."
            },
            "slight right": {
                "default": "약간 우회전하세요.",
                "name": "약간 우회전 하고 {way_name}로 가세요.",
                "destination": "약간 우회전 하신 후 {destination}로 가세요."
            },
            "uturn": {
                "default": "유턴 하세요",
                "name": "유턴해서 {way_name}로 가세요.",
                "destination": "유턴하신 후 {destination}로 가세요."
            }
        },
        "depart": {
            "default": {
                "default": "{direction}로 가세요",
                "name": "{direction} 로 가서 {way_name} 를 이용하세요. ",
                "namedistance": "{direction}로 가서{way_name} 를 {distance}까지 가세요."
            }
        },
        "end of road": {
            "default": {
                "default": "{modifier} 회전하세요.",
                "name": "{modifier}회전하고 {way_name}로 가세요.",
                "destination": "{modifier}회전 하신 후 {destination}로 가세요."
            },
            "straight": {
                "default": "계속 직진해 주세요.",
                "name": "{way_name}로 계속 직진해 주세요.",
                "destination": "{destination}까지 직진해 주세요."
            },
            "uturn": {
                "default": "도로 끝까지 가서 유턴해 주세요.",
                "name": "도로 끝까지 가서 유턴해서 {way_name}로 가세요.",
                "destination": "도로 끝까지 가서 유턴해서 {destination} 까지 가세요."
            }
        },
        "fork": {
            "default": {
                "default": "갈림길에서 {modifier} 으로 가세요.",
                "name": "{modifier}하고 {way_name}로 가세요.",
                "destination": "{modifier}하고 {destination}까지 가세요."
            },
            "slight left": {
                "default": "갈림길에서 좌회전 하세요.",
                "name": "좌회전 해서 {way_name}로 가세요.",
                "destination": "좌회전 해서 {destination}까지 가세요."
            },
            "slight right": {
                "default": "갈림길에서 우회전 하세요.",
                "name": "우회전 해서 {way_name}로 가세요.",
                "destination": "우회전 해서 {destination}까지 가세요."
            },
            "sharp left": {
                "default": "갈림길에서 급좌회전 하세요.",
                "name": "급좌회전 해서 {way_name}로 가세요.",
                "destination": "급좌회전 해서 {destination}까지 가세요."
            },
            "sharp right": {
                "default": "갈림길에서 급우회전 하세요.",
                "name": "급우회전 해서 {way_name}로 가세요.",
                "destination": "급우회전 해서 {destination}까지 가세요."
            },
            "uturn": {
                "default": "유턴하세요.",
                "name": "유턴해서 {way_name}로 가세요.",
                "destination": "유턴해서 {destination}까지 가세요."
            }
        },
        "merge": {
            "default": {
                "default": "{modifier} 합류",
                "name": "{modifier} 합류하여 {way_name}로 가세요.",
                "destination": "{modifier} 합류하여 {destination}로 가세요."
            },
            "straight": {
                "default": "합류",
                "name": "{way_name}로 합류하세요.",
                "destination": "{destination}로 합류하세요."
            },
            "slight left": {
                "default": "좌측으로 합류하세요.",
                "name": "좌측{way_name}로 합류하세요.",
                "destination": "좌측으로 합류하여 {destination}까지 가세요."
            },
            "slight right": {
                "default": "우측으로 합류하세요.",
                "name": "우측{way_name}로 합류하세요.",
                "destination": "우측으로 합류하여 {destination}까지 가세요."
            },
            "sharp left": {
                "default": "좌측으로 합류하세요.",
                "name": "좌측{way_name}로 합류하세요.",
                "destination": "좌측으로 합류하여 {destination}까지 가세요."
            },
            "sharp right": {
                "default": "우측으로 합류하세요.",
                "name": "우측{way_name}로 합류하세요.",
                "destination": "우측으로 합류하여 {destination}까지 가세요."
            },
            "uturn": {
                "default": "유턴하세요.",
                "name": "유턴해서 {way_name}로 가세요.",
                "destination": "유턴해서 {destination}까지 가세요."
            }
        },
        "new name": {
            "default": {
                "default": "{modifier} 유지하세요.",
                "name": "{modifier} 유지해서 {way_name}로 가세요.",
                "destination": "{modifier} 유지해서 {destination}까지 가세요."
            },
            "straight": {
                "default": "직진해주세요.",
                "name": "{way_name}로 계속 가세요.",
                "destination": "{destination}까지 계속 가세요."
            },
            "sharp left": {
                "default": "급좌회전 하세요.",
                "name": "급좌회전 해서 {way_name}로 가세요.",
                "destination": "급좌회전 해서 {destination}까지 가세요."
            },
            "sharp right": {
                "default": "급우회전 하세요.",
                "name": "급우회전 해서 {way_name}로 가세요.",
                "destination": "급우회전 해서 {destination}까지 가세요."
            },
            "slight left": {
                "default": "약간 좌회전 해세요.",
                "name": "약간 좌회전해서 {way_name}로 가세요.",
                "destination": "약간 좌회전 해서 {destination}까지 가세요."
            },
            "slight right": {
                "default": "약간 우회전 해세요.",
                "name": "약간 우회전해서 {way_name}로 가세요.",
                "destination": "약간 우회전 해서 {destination}까지 가세요."
            },
            "uturn": {
                "default": "유턴해주세요.",
                "name": "유턴해서 {way_name}로 가세요.",
                "destination": "유턴해서 {destination}까지 가세요."
            }
        },
        "notification": {
            "default": {
                "default": "{modifier} 하세요.",
                "name": "{modifier}해서 {way_name}로 가세요.",
                "destination": "{modifier}해서 {destination}까지 가세요."
            },
            "uturn": {
                "default": "유턴하세요.",
                "name": "유턴해서 {way_name}로 가세요.",
                "destination": "유턴해서 {destination}까지 가세요."
            }
        },
        "off ramp": {
            "default": {
                "default": "램프로 진출해 주세요..",
                "name": "램프로 진출해서 {way_name}로 가세요.",
                "destination": "램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 출구로 나가세요.",
                "exit_destination": "{exit} 출구로 나가서 {destination}까지 가세요."
            },
            "left": {
                "default": "왼쪽의 램프로 진출해 주세요.",
                "name": "왼쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 왼쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 왼쪽의 출구로 가나서 {destination}까지 가세요."
            },
            "right": {
                "default": "오른쪽의 램프로 진출해 주세요.",
                "name": "오른쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 오른쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 오른쪽의 출구로 가나서 {destination}까지 가세요."
            },
            "sharp left": {
                "default": "왼쪽의 램프로 진출해 주세요.",
                "name": "왼쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 왼쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 왼쪽의 출구로 가나서 {destination}까지 가세요."
            },
            "sharp right": {
                "default": "오른쪽의 램프로 진출해 주세요.",
                "name": "오른쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 오른쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 오른쪽의 출구로 가나서 {destination}까지 가세요."
            },
            "slight left": {
                "default": "왼쪽의 램프로 진출해 주세요.",
                "name": "왼쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 왼쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 왼쪽의 출구로 가나서 {destination}까지 가세요."
            },
            "slight right": {
                "default": "오른쪽의 램프로 진출해 주세요.",
                "name": "오른쪽의 램프로 진출해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진출해서 {destination}까지 가세요.",
                "exit": "{exit} 오른쪽의 출구로 나가세요.",
                "exit_destination": "{exit} 오른쪽의 출구로 가나서 {destination}까지 가세요."
            }
        },
        "on ramp": {
            "default": {
                "default": "램프로 진입해 주세요..",
                "name": "램프로 진입해서 {way_name}로 가세요.",
                "destination": "램프로 진입해서 {destination}까지 가세요."
            },
            "left": {
                "default": "왼쪽의 램프로 진입해 주세요.",
                "name": "왼쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진입해서 {destination}까지 가세요."
            },
            "right": {
                "default": "오른쪽의 램프로 진입해 주세요.",
                "name": "오른쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진입해서 {destination}까지 가세요."
            },
            "sharp left": {
                "default": "왼쪽의 램프로 진입해 주세요.",
                "name": "왼쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진입해서 {destination}까지 가세요."
            },
            "sharp right": {
                "default": "오른쪽의 램프로 진입해 주세요.",
                "name": "오른쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진입해서 {destination}까지 가세요."
            },
            "slight left": {
                "default": "왼쪽의 램프로 진입해 주세요.",
                "name": "왼쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "왼쪽의 램프로 진입해서 {destination}까지 가세요."
            },
            "slight right": {
                "default": "오른쪽의 램프로 진입해 주세요.",
                "name": "오른쪽의 램프로 진입해서 {way_name}로 가세요.",
                "destination": "오른쪽의 램프로 진입해서 {destination}까지 가세요."
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "로터리로 진입하세요.",
                    "name": "로터리로 진입해서 {way_name} 나가세요.",
                    "destination": "로터리로 진입해서 {destination}로 나가세요."
                },
                "name": {
                    "default": "{rotary_name}로 진입하세요.",
                    "name": "{rotary_name}로 진입해서 {way_name}로 나가세요.",
                    "destination": "{rotary_name}로 진입해서 {destination}로 나가세요."
                },
                "exit": {
                    "default": "로터리로 진입해서 {exit_number} 출구로 나가세요.",
                    "name": "로터리로 진입해서 {exit_number} 출구로 나가 {way_name}로 가세요.",
                    "destination": "로터리로 진입해서 {exit_number} 출구로 나가 {destination}로 가세요."
                },
                "name_exit": {
                    "default": "{rotary_name}로 진입해서 {exit_number}번 출구로 나가세요.",
                    "name": "{rotary_name}로 진입해서 {exit_number}번 출구로 나가 {way_name}로 가세요.",
                    "destination": "{rotary_name}로 진입해서 {exit_number}번 출구로 나가 {destination}로 가세요."
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "로터리로 진입해서 {exit_number}로 나가세요.",
                    "name": "로터리로 진입해서 {exit_number}로 나가서 {way_name}로 가세요.",
                    "destination": "로터리로 진입해서 {exit_number}로 나가서 {destination}로 가세요."
                },
                "default": {
                    "default": "로터리로 진입하세요.",
                    "name": "로터리로 진입해서 {way_name} 나가세요.",
                    "destination": "로터리로 진입해서 {destination}로 나가세요."
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier} 하세요.",
                "name": "{modifier} 하시고 {way_name}로 가세요.",
                "destination": "{modifier} 하시고 {destination}까지 가세요."
            },
            "left": {
                "default": "좌회전 하세요.",
                "name": "좌회전 하시고 {way_name}로 가세요.",
                "destination": "좌회전 하시고 {destination}까지 가세요."
            },
            "right": {
                "default": "우회전 하세요.",
                "name": "우회전 하시고 {way_name}로 가세요.",
                "destination": "우회전 하시고 {destination}까지 가세요."
            },
            "straight": {
                "default": "직진 하세요.",
                "name": "직진하시고 {way_name}로 가세요.",
                "destination": "직진하시고 {destination}까지 가세요."
            }
        },
        "exit roundabout": {
            "default": {
                "default": "로타리에서 진출하세요.",
                "name": "로타리에서 진출해서 {way_name}로 가세요.",
                "destination": "로타리에서 진출해서 {destination}까지 가세요."
            }
        },
        "exit rotary": {
            "default": {
                "default": "로타리에서 진출하세요.",
                "name": "로타리에서 진출해서 {way_name}로 가세요.",
                "destination": "로타리에서 진출해서 {destination}까지 가세요."
            }
        },
        "turn": {
            "default": {
                "default": "{modifier} 하세요.",
                "name": "{modifier} 하시고 {way_name}로 가세요.",
                "destination": "{modifier} 하시고 {destination}까지 가세요."
            },
            "left": {
                "default": "좌회전 하세요.",
                "name": "좌회전 하시고 {way_name}로 가세요.",
                "destination": "좌회전 하시고 {destination}까지 가세요."
            },
            "right": {
                "default": "우회전 하세요.",
                "name": "우회전 하시고 {way_name}로 가세요.",
                "destination": "우회전 하시고 {destination}까지 가세요."
            },
            "straight": {
                "default": "직진 하세요.",
                "name": "직진하시고 {way_name}로 가세요.",
                "destination": "직진하시고 {destination}까지 가세요."
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "직진하세요."
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],35:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": false
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "ပထမ",
                "2": "ဒုတိယ",
                "3": "တတိယ",
                "4": "စတုတၳ",
                "5": "ပဥၥမ",
                "6": "ဆဌမ",
                "7": "သတၱမ",
                "8": "အဌမ",
                "9": "နဝမ",
                "10": "ဒသမ"
            },
            "direction": {
                "north": "ေျမာက္အရပ္",
                "northeast": "အေရွ႕ေျမာက္အရပ္",
                "east": "အေရွ႕အရပ္",
                "southeast": "အေရွ႕ေတာင္အရပ္",
                "south": "ေတာင္အရပ္",
                "southwest": "အေနာက္ေတာင္အရပ္",
                "west": "အေနာက္အရပ္",
                "northwest": "အေနာက္ေျမာက္အရပ္"
            },
            "modifier": {
                "left": "ဘယ္ဘက္",
                "right": "ညာဘက္",
                "sharp left": "ဘယ္ဘက္ ေထာင့္ခ်ိဳး",
                "sharp right": "ညာဘက္ ေထာင္႔ခ်ိဳး",
                "slight left": "ဘယ္ဘက္ အနည္းငယ္",
                "slight right": "ညာဘက္ အနည္းငယ္",
                "straight": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း",
                "uturn": "ဂ-ေကြ႔"
            },
            "lanes": {
                "xo": "ညာဘက္သို႕ဆက္သြားပါ",
                "ox": "ဘယ္ဘက္သို႕ဆက္သြားပါ",
                "xox": "အလယ္တြင္ဆက္ေနပါ",
                "oxo": "ဘယ္ သို႕မဟုတ္ ညာဘက္သို႕ ဆက္သြားပါ"
            }
        },
        "modes": {
            "ferry": {
                "default": "ဖယ္ရီ စီးသြားပါ",
                "name": "{way_name}ကို ဖယ္ရီစီးသြားပါ",
                "destination": "{destination}ဆီသို႕ ဖယ္ရီစီးသြားပါ"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}ျပီးေနာက္ {distance}အတြင္း {instruction_two}",
            "two linked": "{instruction_one}ျပီးေနာက္ {instruction_two}",
            "one in distance": "{distance}အတြင္း {instruction_one}",
            "name and ref": "{name}( {ref})",
            "exit with number": "{exit}မွထြက္ပါ"
        },
        "arrive": {
            "default": {
                "default": "{nth}သင္ သြားလိုေသာ ခရီးပန္းတိုင္သို႕ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ေရာက္လိမ့္မည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name} မွာ ေရာက္ရွိျပီ"
            },
            "left": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္လိမ့္မည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name}မွာဘယ္ဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "right": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ{nth} ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္လိမ့္မည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name} မွာညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "sharp left": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္ရွိျပီ",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name}မွာဘယ္ဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "sharp right": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ{nth} ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္လိမ့္မည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name} မွာညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "slight right": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ{nth} ခရီးပန္းတိုင္သို႕ ညာဘက္ေကြ႕ကာ ေရာက္လိမ့္မည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရိွၿပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name} မွာညာဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "slight left": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕ဘယ္ဘက္တြင္ေရာက္ရွိျပီ",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရွိျပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name}မွာဘယ္ဘက္ေကြ႕ကာ ေရာက္ရွိျပီ"
            },
            "straight": {
                "default": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕တည့္တည့္သြားကာရာက္ရွိျပီ",
                "upcoming": "သင္ သြားလိုေသာ {nth}ခရီးပန္းတိုင္သို႕တည့္တည့္သြားကာရာက္ရွိမည္",
                "short": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္ရွိျပီ",
                "short-upcoming": "သင္သြားလိုေသာ ေနရာသို႔ ေရာက္လိမ့္မည္",
                "named": "သင္ သည္ {waypoint_name}မွာတည့္တည့္သြားကာ ေရာက္ရွိျပီ"
            }
        },
        "continue": {
            "default": {
                "default": "{modifier}ကိုလွည့္ပါ",
                "name": "{way_name}​​ေပၚတြင္ေနရန္ {modifier}ကိုလွည့္ပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို လွည္႕ပါ",
                "exit": "{way_name}​​ေပၚသို႕ {modifier}ကိုလွည့္ပါ"
            },
            "straight": {
                "default": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚတြင္ေနရန္တည္တည့္ဆက္သြာပါ",
                "destination": "{destination}ဆီသို႕ဆက္သြားပါ",
                "distance": "{distance}ေလာက္ တည့္တည့္ ဆက္သြားပါ",
                "namedistance": "{way_name}​​ေပၚတြင္{distance}ေလာက္ဆက္သြားပါ"
            },
            "sharp left": {
                "default": "ဘယ္ဘက္ေထာင့္ခ်ိဳးေကြ႕ပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ဘယ္ဘက္ေထာင့္ခ်ိဳးေကြ႕ပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္ေထာင့္ခ်ိဳးေကြ႕ပါ"
            },
            "sharp right": {
                "default": "ညာဘက္ ေထာင္႔ခ်ိဳးေကြ႕ပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ညာဘက္ေထာင့္ခ်ိဳးေကြ႕ပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္ေထာင့္ခ်ိဳးေကြ႕ပါ"
            },
            "slight left": {
                "default": "ဘယ္ဘက္ အနည္းငယ္ေကြ႕ပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ဘယ္ဘက္အနည္းငယ္ေကြ႕ပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္အနည္းငယ္ခ်ိဳးေကြ႕ပါ"
            },
            "slight right": {
                "default": "ညာဘက္ အနည္းငယ္ခ်ိဳးေကြ႕ပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ညာဘက္အနည္းငယ္ေကြ႕ပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္အနည္းငယ္ခ်ိဳးေကြ႕ပါ"
            },
            "uturn": {
                "default": "ဂ-ေကြ႔ ေကြ႔ပါ",
                "name": "{way_name}လမ္းဘက္သို႕ ဂ-ေကြ႕ေကြ႕ျပီးဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "depart": {
            "default": {
                "default": "{direction}သို႕ ဦးတည္ပါ",
                "name": "{direction}ကို {way_name}အေပၚတြင္ ဦးတည္ပါ",
                "namedistance": "{direction}ကို {way_name}အေပၚတြင္{distance}ေလာက္ ဦးတည္ဆက္သြားပါ"
            }
        },
        "end of road": {
            "default": {
                "default": "{modifier}သို႕လွည့္ပါ",
                "name": "{way_name}​​ေပၚသို႕ {modifier}ကိုလွည့္ပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို လွည္႕ပါ"
            },
            "straight": {
                "default": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕တည့္တည့္ဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕တည့္တည့္ဆက္သြားပါ"
            },
            "uturn": {
                "default": "လမ္းအဆံုးတြင္ ဂ-ေကြ႕ေကြ႕ပါ",
                "name": "လမ္းအဆံုးတြင္ {way_name}​​ေပၚသို႕ဂ-ေကြ႕ေကြ႕ပါ",
                "destination": "လမ္းအဆံုးတြင္{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "fork": {
            "default": {
                "default": "လမ္းဆံုလမ္းခြတြင္ {modifier}ကိုဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ {modifier}ကိုဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို ဆက္သြားပါ"
            },
            "slight left": {
                "default": "လမ္းဆံုလမ္းခြတြင္ဘယ္ဘက္ကိုဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ဘယ္ဘက္ကိုဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္ကို ဆက္သြားပါ"
            },
            "slight right": {
                "default": "လမ္းဆံုလမ္းခြတြင္ညာဘက္ကိုဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ညာဘက္ကိုဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္ကို ဆက္သြားပါ"
            },
            "sharp left": {
                "default": "လမ္းဆံုလမ္းခြတြင္ဘယ္ဘက္ေထာင့္ခ်ိဳးကိုသြားပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ဘယ္ဘက္ေထာင့္ခ်ိဳးယူပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္ေထာင့္ခ်ိဳး သြားပါ"
            },
            "sharp right": {
                "default": "လမ္းဆံုလမ္းခြတြင္ညာဘက္ေထာင့္ခ်ိဳးကိုသြားပါ",
                "name": "{way_name}​ေပၚသို႕ ညာဘက္ေထာင့္ခ်ိဳးယူပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္ေထာင့္ခ်ိဳး သြားပါ"
            },
            "uturn": {
                "default": "ဂ-ေကြ႔ ေကြ႔ပါ",
                "name": "{way_name}သို႕ဂ-ေကြ႕ေကြ႕ပါ",
                "destination": "{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "merge": {
            "default": {
                "default": "{modifier}ကိုလာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕ {modifier}ကိုလာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို လာေရာက္ေပါင္းဆံုပါ"
            },
            "straight": {
                "default": "လာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕လာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ လာေရာက္ေပါင္းဆံုပါ"
            },
            "slight left": {
                "default": "ဘယ္ဘက္သို႕လာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕ဘယ္ဘက္ကိုလာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္ကို လာေရာက္ေပါင္းဆံုပါ"
            },
            "slight right": {
                "default": "ညာဘက္သို႕လာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕ညာဘက္ကိုလာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္ကို လာေရာက္ေပါင္းဆံုပါ"
            },
            "sharp left": {
                "default": "ဘယ္ဘက္သို႕လာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕ဘယ္ဘက္ကိုလာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္ကို လာေရာက္ေပါင္းဆံုပါ"
            },
            "sharp right": {
                "default": "ညာဘက္သို႕လာေရာက္ေပါင္းဆံုပါ",
                "name": "{way_name}​​ေပၚသို႕ညာဘက္ကိုလာေရာက္ေပါင္းဆံုပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္ကို လာေရာက္ေပါင္းဆံုပါ"
            },
            "uturn": {
                "default": "ဂ-ေကြ႔ ေကြ႕ပါ",
                "name": "{way_name}လမ္းဘက္သို႔ ဂ-ေကြ႔ ေကြ႔ပါ ",
                "destination": "{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "new name": {
            "default": {
                "default": "{modifier}ကိုဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ {modifier}ကိုဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို ဆက္သြားပါ"
            },
            "straight": {
                "default": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ဆက္သြားပါ"
            },
            "sharp left": {
                "default": "ဘယ္ဘက္ေထာင့္ခ်ိဳးယူပါ",
                "name": "{way_name}​ေပၚတြင္ေနရန္ ဘယ္ဘက္ေထာင့္ခ်ိဳးယူပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္ေထာင့္ခ်ိဳး သြားပါ"
            },
            "sharp right": {
                "default": "ညာဘက္ ေထာင္႔ခ်ိဳးယူပါ",
                "name": "{way_name}​ေပၚသို႕ ညာဘက္ေထာင့္ခ်ိဳးယူပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္ေထာင့္ခ်ိဳး သြားပါ"
            },
            "slight left": {
                "default": "ဘယ္ဘက္ အနည္းငယ္ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ဘယ္ဘက္ အနည္းငယ္ဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္အနည္းငယ္ဆက္သြားပါ"
            },
            "slight right": {
                "default": "ညာဘက္ အနည္းငယ္ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ညာဘက္ အနည္းငယ္ဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ညာဘက္အနည္းငယ္ဆက္သြားပါ"
            },
            "uturn": {
                "default": "ဂ-ေကြ႔ ေကြ႔ပါ",
                "name": "{way_name}လမ္းဘက္သို႔ ဂ-ေကြ႔ ေကြ႔ပါ",
                "destination": "{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "notification": {
            "default": {
                "default": "{modifier}ကိုဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕ {modifier}ကိုဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕ {modifier}ကို ဆက္သြားပါ"
            },
            "uturn": {
                "default": "ဂ-ေကြ႔ ေကြ႔ပါ",
                "name": "{way_name}လမ္းဘက္သို႔ ဂ-ေကြ႔ ေကြ႔ပါ",
                "destination": "{destination}ဆီသို႕ ဂေကြ႕ခ်ိဳးေကြ႕ပါ"
            }
        },
        "off ramp": {
            "default": {
                "default": "ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ {exit} ကိုယူပါ"
            },
            "left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ဘယ္ဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ဘယ္ဘက္မွ {exit} ကိုယူပါ"
            },
            "right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ညာဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ညာဘက္မွ {exit} ကိုယူပါ"
            },
            "sharp left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ဘယ္ဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ဘယ္ဘက္မွ {exit} ကိုယူပါ"
            },
            "sharp right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ညာဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ညာဘက္မွ {exit} ကိုယူပါ"
            },
            "slight left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ဘယ္ဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ဘယ္ဘက္မွ {exit} ကိုယူပါ"
            },
            "slight right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "exit": "ညာဘက္တြင္{exit}ကို ယူပါ",
                "exit_destination": "{destination}ဆီသို႕ညာဘက္မွ {exit} ကိုယူပါ"
            }
        },
        "on ramp": {
            "default": {
                "default": "ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "sharp left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "sharp right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "slight left": {
                "default": "ဘယ္ဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ဘယ္ဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            },
            "slight right": {
                "default": "ညာဘက္သို႕ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ ​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ",
                "destination": "{destination}ဆီသို႕ ညာဘက္​ေပၚတြင္ခ်ဥ္းကပ္လမ္းကိုယူပါ"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "အဝိုင္းပတ္သို႕ဝင္ပါ",
                    "name": "{way_name}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ ",
                    "destination": "{destination}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ"
                },
                "name": {
                    "default": "{rotary_name}သို႕ဝင္ပါ",
                    "name": "{rotary_name}အဝိုင္းပတ္ဝင္ျပီး{way_name}ေပၚသို႕ထြက္ပါ",
                    "destination": "{rotary_name}အဝိုင္းပတ္ဝင္ျပီး{destination}ဆီသို႕ထြက္ပါ"
                },
                "exit": {
                    "default": "အဝိုင္းပတ္ဝင္ျပီး{exit_number}ကိုယူကာျပန္ထြက္ပါ",
                    "name": "အဝိုင္းပတ္သို႕ဝင္ျပီး{exit_number}ကိုယူကာ{way_name}ေပၚသို႕ထြက္ပါ",
                    "destination": "အဝိုင္းပတ္ဝင္ျပီး{exit_number}ကိုယူကာ{destination}ဆီသို႕ထြက္ပါ"
                },
                "name_exit": {
                    "default": "{rotary_name}ကိုဝင္ျပီး {exit_number}ကိုယူကာထြက္ပါ",
                    "name": "{rotary_name}ကိုဝင္ျပီး{exit_number}ကိုယူကာ{way_name}ေပၚသို႕ထြက္ပါ",
                    "destination": "{rotary_name}ဝင္ျပီး{exit_number}ကိုယူကာ{destination}ဆီသို႕ထြက္ပါ"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "{exit_number}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                    "name": "အဝိုင္းပတ္ဝင္ျပီး{exit_number}ကိုယူကာ{way_name}ေပၚသို႕ထြက္ပါ",
                    "destination": "အဝိုင္းပတ္ဝင္ျပီး{exit_number}ကိုယူကာ{destination}ဆီသို႕ထြက္ပါ"
                },
                "default": {
                    "default": "အဝိုင္းပတ္ဝင္ပါ",
                    "name": "{way_name}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                    "destination": "{destination}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier}ကိုလွည့္ပါ ",
                "name": "{modifier}​ေပၚသို{way_name}ကိုဆက္သြားပါ ",
                "destination": "{modifier}ဆီသို႕{destination}ကို ဆက္သြားပါ "
            },
            "left": {
                "default": "ဘယ္ဘက္သို႕ျပန္လွည္႔ပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ကိုဆက္သြားပါ ",
                "destination": "{destination}ဆီသို႕ဘယ္ဘက္မွ ေကြ႔ပါ"
            },
            "right": {
                "default": "ညာဘက္သို႔ျပန္လွည္႔ပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ကိုလာေရာက္ေပါင္းဆံုပါ ",
                "destination": "{destination}ညာဘက္သို႔ ေကြ႔ပါ"
            },
            "straight": {
                "default": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း ဆက္သြားပါ",
                "name": "{way_name}​​ေပၚသို႕တည့္တည့္ဆက္သြားပါ",
                "destination": "{destination}ဆီသို႕တည့္တည့္ဆက္သြားပါ"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                "name": "{way_name}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                "destination": "ဦးတည္အဝိုင္းပတ္လမ္းမွထြက္ပါ{destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "အဝိုင္းပတ္လမ္းမွထြက္ပါဦးတည္အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                "name": "{way_name}ေပၚသို႔အဝိုင္းပတ္လမ္းမွထြက္ပါ",
                "destination": "ဦးတည္အဝိုင္းပတ္လမ္းမွထြက္ပါ{destination}"
            }
        },
        "turn": {
            "default": {
                "default": "{modifier}ကိုလွည့္ပါ ",
                "name": "{modifier}​ေပၚသို{way_name}ကိုဆက္သြားပါ ",
                "destination": "{modifier}ဆီသို႕{destination}ကို ဆက္သြားပါ "
            },
            "left": {
                "default": "ဘယ္ဘက္သို႕ျပန္လွည္႔ပါ",
                "name": "{way_name}​ေပၚသို႕ဘယ္ဘက္ကိုဆက္သြားပါ ",
                "destination": "{destination}ဘယ္ဘက္သို႔ ေကြ႔ပါ"
            },
            "right": {
                "default": "ညာဘက္သို႔ျပန္လွည္႔ပါ",
                "name": "{way_name}​ေပၚသို႕ညာဘက္ကိုလာေရာက္ေပါင္းဆံုပါ ",
                "destination": "{destination}ညာဘက္သို႔ ေကြ႔ပါ"
            },
            "straight": {
                "default": "တည္႔တည္႔သြားပါ",
                "name": "{way_name}",
                "destination": "{destination}ဆီသို႕တည့္တည့္သြားပါ"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "ေျဖာင္႔ေျဖာင္႔တန္းတန္း ဆက္သြားပါ"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],36:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1e",
                "2": "2e",
                "3": "3e",
                "4": "4e",
                "5": "5e",
                "6": "6e",
                "7": "7e",
                "8": "8e",
                "9": "9e",
                "10": "10e"
            },
            "direction": {
                "north": "noord",
                "northeast": "noordoost",
                "east": "oost",
                "southeast": "zuidoost",
                "south": "zuid",
                "southwest": "zuidwest",
                "west": "west",
                "northwest": "noordwest"
            },
            "modifier": {
                "left": "links",
                "right": "rechts",
                "sharp left": "scherpe bocht naar links",
                "sharp right": "scherpe bocht naar rechts",
                "slight left": "iets naar links",
                "slight right": "iets naar rechts",
                "straight": "rechtdoor",
                "uturn": "omkeren"
            },
            "lanes": {
                "xo": "Rechts aanhouden",
                "ox": "Links aanhouden",
                "xox": "In het midden blijven",
                "oxo": "Links of rechts blijven"
            }
        },
        "modes": {
            "ferry": {
                "default": "Neem de veerpont",
                "name": "Neem de veerpont {way_name}",
                "destination": "Neem de veerpont richting {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, dan na {distance}, {instruction_two}",
            "two linked": "{instruction_one}, daarna {instruction_two}",
            "one in distance": "Over {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "afslag {exit}"
        },
        "arrive": {
            "default": {
                "default": "Je bent gearriveerd op de {nth} bestemming.",
                "upcoming": "U arriveert op de {nth} bestemming",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name}"
            },
            "left": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich links.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de linkerkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name}, de bestemming is aan de linkerkant"
            },
            "right": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich rechts.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de rechterkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name}, de bestemming is aan de  rechterkant"
            },
            "sharp left": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich links.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de linkerkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name}, de bestemming is aan de linkerkant"
            },
            "sharp right": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich rechts.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de rechterkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name},  de bestemming is aan de rechterkant"
            },
            "slight right": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich rechts.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de rechterkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name},  de bestemming is aan de rechterkant"
            },
            "slight left": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich links.",
                "upcoming": "Uw {nth} bestemming bevindt zich aan de linkerkant",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name},  de bestemming is aan de linkerkant"
            },
            "straight": {
                "default": "Je bent gearriveerd. De {nth} bestemming bevindt zich voor je.",
                "upcoming": "Uw {nth} bestemming is recht voor u",
                "short": "U bent gearriveerd",
                "short-upcoming": "U zult aankomen",
                "named": "U bent gearriveerd bij {waypoint_name}, de bestemming is recht voor u"
            }
        },
        "continue": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Sla {modifier} om op {way_name} te blijven",
                "destination": "Ga {modifier} richting {destination}",
                "exit": "Ga {modifier} naar {way_name}"
            },
            "straight": {
                "default": "Ga rechtdoor",
                "name": "Blijf rechtdoor gaan op {way_name}",
                "destination": "Ga rechtdoor richting {destination}",
                "distance": "Ga rechtdoor voor {distance}",
                "namedistance": "Ga verder op {way_name} voor {distance}"
            },
            "sharp left": {
                "default": "Linksaf",
                "name": "Sla scherp links af om op {way_name} te blijven",
                "destination": "Linksaf richting {destination}"
            },
            "sharp right": {
                "default": "Rechtsaf",
                "name": "Sla scherp rechts af om op {way_name} te blijven",
                "destination": "Rechtsaf richting {destination}"
            },
            "slight left": {
                "default": "Ga links",
                "name": "Links afbuigen om op {way_name} te blijven",
                "destination": "Rechts afbuigen om op {destination} te blijven"
            },
            "slight right": {
                "default": "Rechts afbuigen",
                "name": "Rechts afbuigen om op {way_name} te blijven",
                "destination": "Rechts afbuigen richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Draai om en ga verder op {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Vertrek in {direction}elijke richting",
                "name": "Neem {way_name} in {direction}elijke richting",
                "namedistance": "Ga richting {direction} op {way_name} voor {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Ga {modifier} naar {way_name}",
                "destination": "Ga {modifier} richting {destination}"
            },
            "straight": {
                "default": "Ga in de aangegeven richting",
                "name": "Ga naar {way_name}",
                "destination": "Ga richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Keer om naar {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "Ga {modifier} op de splitsing",
                "name": "Houd {modifier} aan, tot {way_name}",
                "destination": "Houd {modifier}, in de richting van {destination}"
            },
            "slight left": {
                "default": "Links aanhouden op de splitsing",
                "name": "Houd links aan, tot {way_name}",
                "destination": "Houd links aan, richting {destination}"
            },
            "slight right": {
                "default": "Rechts aanhouden op de splitsing",
                "name": "Houd rechts aan, tot {way_name}",
                "destination": "Houd rechts aan, richting {destination}"
            },
            "sharp left": {
                "default": "Neem bij de splitsing, een scherpe bocht, naar links ",
                "name": "Neem een scherpe bocht naar links, tot aan {way_name}",
                "destination": "Neem een scherpe bocht naar links, richting {destination}"
            },
            "sharp right": {
                "default": "Neem  op de splitsing, een scherpe bocht, naar rechts",
                "name": "Neem een scherpe bocht naar rechts, tot aan {way_name}",
                "destination": "Neem een scherpe bocht naar rechts, richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Keer om naar {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Bij de splitsing {modifier}",
                "name": "Bij de splitsing {modifier} naar {way_name}",
                "destination": "Bij de splitsing {modifier} richting {destination}"
            },
            "straight": {
                "default": "Samenvoegen",
                "name": "Ga verder op {way_name}",
                "destination": "Ga verder richting {destination}"
            },
            "slight left": {
                "default": "Bij de splitsing links aanhouden",
                "name": "Bij de splitsing links aanhouden naar {way_name}",
                "destination": "Bij de splitsing links aanhouden richting {destination}"
            },
            "slight right": {
                "default": "Bij de splitsing rechts aanhouden",
                "name": "Bij de splitsing rechts aanhouden naar {way_name}",
                "destination": "Bij de splitsing rechts aanhouden richting {destination}"
            },
            "sharp left": {
                "default": "Bij de splitsing linksaf",
                "name": "Bij de splitsing linksaf naar {way_name}",
                "destination": "Bij de splitsing linksaf richting {destination}"
            },
            "sharp right": {
                "default": "Bij de splitsing rechtsaf",
                "name": "Bij de splitsing rechtsaf naar {way_name}",
                "destination": "Bij de splitsing rechtsaf richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Keer om naar {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Ga {modifier} naar {way_name}",
                "destination": "Ga {modifier} richting {destination}"
            },
            "straight": {
                "default": "Ga in de aangegeven richting",
                "name": "Ga rechtdoor naar {way_name}",
                "destination": "Ga rechtdoor richting {destination}"
            },
            "sharp left": {
                "default": "Neem een scherpe bocht, naar links",
                "name": "Linksaf naar {way_name}",
                "destination": "Linksaf richting {destination}"
            },
            "sharp right": {
                "default": "Neem een scherpe bocht, naar rechts",
                "name": "Rechtsaf naar {way_name}",
                "destination": "Rechtsaf richting {destination}"
            },
            "slight left": {
                "default": "Links aanhouden",
                "name": "Links aanhouden naar {way_name}",
                "destination": "Links aanhouden richting {destination}"
            },
            "slight right": {
                "default": "Rechts aanhouden",
                "name": "Rechts aanhouden naar {way_name}",
                "destination": "Rechts aanhouden richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Keer om naar {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Ga {modifier} naar {way_name}",
                "destination": "Ga {modifier} richting {destination}"
            },
            "uturn": {
                "default": "Keer om",
                "name": "Keer om naar {way_name}",
                "destination": "Keer om richting {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Neem de afrit",
                "name": "Neem de afrit naar {way_name}",
                "destination": "Neem de afrit richting {destination}",
                "exit": "Neem afslag {exit}",
                "exit_destination": "Neem afslag {exit} richting {destination}"
            },
            "left": {
                "default": "Neem de afrit links",
                "name": "Neem de afrit links naar {way_name}",
                "destination": "Neem de afrit links richting {destination}",
                "exit": "Neem afslag {exit} aan de linkerkant",
                "exit_destination": "Neem afslag {exit} aan de linkerkant richting {destination}"
            },
            "right": {
                "default": "Neem de afrit rechts",
                "name": "Neem de afrit rechts naar {way_name}",
                "destination": "Neem de afrit rechts richting {destination}",
                "exit": "Neem afslag {exit} aan de rechterkant",
                "exit_destination": "Neem afslag {exit} aan de rechterkant richting {destination}"
            },
            "sharp left": {
                "default": "Neem de afrit links",
                "name": "Neem de afrit links naar {way_name}",
                "destination": "Neem de afrit links richting {destination}",
                "exit": "Neem afslag {exit} aan de linkerkant",
                "exit_destination": "Neem afslag {exit} aan de linkerkant richting {destination}"
            },
            "sharp right": {
                "default": "Neem de afrit rechts",
                "name": "Neem de afrit rechts naar {way_name}",
                "destination": "Neem de afrit rechts richting {destination}",
                "exit": "Neem afslag {exit} aan de rechterkant",
                "exit_destination": "Neem afslag {exit} aan de rechterkant richting {destination}"
            },
            "slight left": {
                "default": "Neem de afrit links",
                "name": "Neem de afrit links naar {way_name}",
                "destination": "Neem de afrit links richting {destination}",
                "exit": "Neem afslag {exit} aan de linkerkant",
                "exit_destination": "Neem afslag {exit} aan de linkerkant richting {destination}"
            },
            "slight right": {
                "default": "Neem de afrit rechts",
                "name": "Neem de afrit rechts naar {way_name}",
                "destination": "Neem de afrit rechts richting {destination}",
                "exit": "Neem afslag {exit} aan de rechterkant",
                "exit_destination": "Neem afslag {exit} aan de rechterkant richting {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Neem de oprit",
                "name": "Neem de oprit naar {way_name}",
                "destination": "Neem de oprit richting {destination}"
            },
            "left": {
                "default": "Neem de oprit links",
                "name": "Neem de oprit links naar {way_name}",
                "destination": "Neem de oprit links richting {destination}"
            },
            "right": {
                "default": "Neem de oprit rechts",
                "name": "Neem de oprit rechts naar {way_name}",
                "destination": "Neem de oprit rechts richting {destination}"
            },
            "sharp left": {
                "default": "Neem de oprit links",
                "name": "Neem de oprit links naar {way_name}",
                "destination": "Neem de oprit links richting {destination}"
            },
            "sharp right": {
                "default": "Neem de oprit rechts",
                "name": "Neem de oprit rechts naar {way_name}",
                "destination": "Neem de oprit rechts richting {destination}"
            },
            "slight left": {
                "default": "Neem de oprit links",
                "name": "Neem de oprit links naar {way_name}",
                "destination": "Neem de oprit links richting {destination}"
            },
            "slight right": {
                "default": "Neem de oprit rechts",
                "name": "Neem de oprit rechts naar {way_name}",
                "destination": "Neem de oprit rechts richting {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Betreedt de rotonde",
                    "name": "Betreedt rotonde en sla af op {way_name}",
                    "destination": "Betreedt rotonde en sla af richting {destination}"
                },
                "name": {
                    "default": "Ga het knooppunt {rotary_name} op",
                    "name": "Verlaat het knooppunt {rotary_name} naar {way_name}",
                    "destination": "Verlaat het knooppunt {rotary_name} richting {destination}"
                },
                "exit": {
                    "default": "Betreedt rotonde en neem afslag {exit_number}",
                    "name": "Betreedt rotonde en neem afslag {exit_number} naar {way_name}",
                    "destination": "Betreedt rotonde en neem afslag {exit_number} richting {destination}"
                },
                "name_exit": {
                    "default": "Ga het knooppunt {rotary_name} op en neem afslag {exit_number}",
                    "name": "Ga het knooppunt {rotary_name} op en neem afslag {exit_number} naar {way_name}",
                    "destination": "Ga het knooppunt {rotary_name} op en neem afslag {exit_number} richting {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Betreedt rotonde en neem afslag {exit_number}",
                    "name": "Betreedt rotonde en neem afslag {exit_number} naar {way_name}",
                    "destination": "Betreedt rotonde en neem afslag {exit_number} richting {destination}"
                },
                "default": {
                    "default": "Betreedt de rotonde",
                    "name": "Betreedt rotonde en sla af op {way_name}",
                    "destination": "Betreedt rotonde en sla af richting {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Ga {modifier} naar {way_name}",
                "destination": "Ga {modifier} richting {destination}"
            },
            "left": {
                "default": "Ga linksaf",
                "name": "Ga linksaf naar {way_name}",
                "destination": "Ga linksaf richting {destination}"
            },
            "right": {
                "default": "Ga rechtsaf",
                "name": "Ga rechtsaf naar {way_name}",
                "destination": "Ga rechtsaf richting {destination}"
            },
            "straight": {
                "default": "Ga in de aangegeven richting",
                "name": "Ga naar {way_name}",
                "destination": "Ga richting {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Verlaat de rotonde",
                "name": "Verlaat de rotonde en ga verder op {way_name}",
                "destination": "Verlaat de rotonde richting {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Verlaat de rotonde",
                "name": "Verlaat de rotonde en ga verder op {way_name}",
                "destination": "Verlaat de rotonde richting {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Ga {modifier}",
                "name": "Ga {modifier} naar {way_name}",
                "destination": "Ga {modifier} richting {destination}"
            },
            "left": {
                "default": "Ga linksaf",
                "name": "Ga linksaf naar {way_name}",
                "destination": "Ga linksaf richting {destination}"
            },
            "right": {
                "default": "Ga rechtsaf",
                "name": "Ga rechtsaf naar {way_name}",
                "destination": "Ga rechtsaf richting {destination}"
            },
            "straight": {
                "default": "Ga rechtdoor",
                "name": "Ga rechtdoor naar {way_name}",
                "destination": "Ga rechtdoor richting {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Rechtdoor"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],37:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1.",
                "2": "2.",
                "3": "3.",
                "4": "4.",
                "5": "5.",
                "6": "6.",
                "7": "7.",
                "8": "8.",
                "9": "9.",
                "10": "10."
            },
            "direction": {
                "north": "nord",
                "northeast": "nordøst",
                "east": "øst",
                "southeast": "sørøst",
                "south": "sør",
                "southwest": "sørvest",
                "west": "vest",
                "northwest": "nordvest"
            },
            "modifier": {
                "left": "venstre",
                "right": "høyre",
                "sharp left": "skarp venstre",
                "sharp right": "skarp høyre",
                "slight left": "litt til venstre",
                "slight right": "litt til høyre",
                "straight": "rett frem",
                "uturn": "U-sving"
            },
            "lanes": {
                "xo": "Hold til høyre",
                "ox": "Hold til venstre",
                "xox": "Hold deg i midten",
                "oxo": "Hold til venstre eller høyre"
            }
        },
        "modes": {
            "ferry": {
                "default": "Ta ferja",
                "name": "Ta ferja {way_name}",
                "destination": "Ta ferja til {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, deretter {instruction_two} om {distance}",
            "two linked": "{instruction_one}, deretter {instruction_two}",
            "one in distance": "Om {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "avkjørsel {exit}"
        },
        "arrive": {
            "default": {
                "default": "Du har ankommet din {nth} destinasjon",
                "upcoming": "Du vil ankomme din {nth} destinasjon",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}"
            },
            "left": {
                "default": "Du har ankommet din {nth} destinasjon, på din venstre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din venstre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din venstre side"
            },
            "right": {
                "default": "Du har ankommet din {nth} destinasjon, på din høyre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din høyre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din høyre side"
            },
            "sharp left": {
                "default": "Du har ankommet din {nth} destinasjon, på din venstre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din venstre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din venstre side"
            },
            "sharp right": {
                "default": "Du har ankommet din {nth} destinasjon, på din høyre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din høyre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din høyre side"
            },
            "slight right": {
                "default": "Du har ankommet din {nth} destinasjon, på din høyre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din høyre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din høyre side"
            },
            "slight left": {
                "default": "Du har ankommet din {nth} destinasjon, på din venstre side",
                "upcoming": "Du vil ankomme din {nth} destinasjon, på din venstre side",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, på din venstre side"
            },
            "straight": {
                "default": "Du har ankommet din {nth} destinasjon, rett forut",
                "upcoming": "Du vil ankomme din {nth} destinasjon, rett forut",
                "short": "Du har ankommet",
                "short-upcoming": "Du vil ankomme",
                "named": "Du har ankommet {waypoint_name}, rett forut"
            }
        },
        "continue": {
            "default": {
                "default": "Ta til {modifier}",
                "name": "Ta til {modifier} for å bli værende på {way_name}",
                "destination": "Ta til {modifier} mot {destination}",
                "exit": "Ta til {modifier} inn på {way_name}"
            },
            "straight": {
                "default": "Fortsett rett frem",
                "name": "Fortsett rett frem for å bli værende på {way_name}",
                "destination": "Fortsett mot {destination}",
                "distance": "Fortsett rett frem, {distance} ",
                "namedistance": "Fortsett på {way_name}, {distance}"
            },
            "sharp left": {
                "default": "Sving skarpt til venstre",
                "name": "Sving skarpt til venstre for å bli værende på {way_name}",
                "destination": "Sving skarpt til venstre mot {destination}"
            },
            "sharp right": {
                "default": "Sving skarpt til høyre",
                "name": "Sving skarpt til høyre for å bli værende på {way_name}",
                "destination": "Sving skarpt mot {destination}"
            },
            "slight left": {
                "default": "Sving svakt til venstre",
                "name": "Sving svakt til venstre for å bli værende på {way_name}",
                "destination": "Sving svakt til venstre mot {destination}"
            },
            "slight right": {
                "default": "Sving svakt til høyre",
                "name": "Sving svakt til høyre for å bli værende på {way_name}",
                "destination": "Sving svakt til høyre mot {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving",
                "name": "Ta en U-sving og fortsett på {way_name}",
                "destination": "Ta en U-sving mot {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Kjør i retning {direction}",
                "name": "Kjør i retning {direction} på {way_name}",
                "namedistance": "Kjør i retning {direction} på {way_name}, {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Sving {modifier}",
                "name": "Ta til {modifier} inn på {way_name}",
                "destination": "Sving {modifier} mot {destination}"
            },
            "straight": {
                "default": "Fortsett rett frem",
                "name": "Fortsett rett frem til  {way_name}",
                "destination": "Fortsett rett frem mot {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving i enden av veien",
                "name": "Ta en U-sving til {way_name} i enden av veien",
                "destination": "Ta en U-sving mot {destination} i enden av veien"
            }
        },
        "fork": {
            "default": {
                "default": "Hold til {modifier} i veikrysset",
                "name": "Hold til {modifier} inn på {way_name}",
                "destination": "Hold til {modifier} mot {destination}"
            },
            "slight left": {
                "default": "Hold til venstre i veikrysset",
                "name": "Hold til venstre inn på {way_name}",
                "destination": "Hold til venstre mot {destination}"
            },
            "slight right": {
                "default": "Hold til høyre i veikrysset",
                "name": "Hold til høyre inn på {way_name}",
                "destination": "Hold til høyre mot {destination}"
            },
            "sharp left": {
                "default": "Sving skarpt til venstre i veikrysset",
                "name": "Sving skarpt til venstre inn på {way_name}",
                "destination": "Sving skarpt til venstre mot {destination}"
            },
            "sharp right": {
                "default": "Sving skarpt til høyre i veikrysset",
                "name": "Sving skarpt til høyre inn på {way_name}",
                "destination": "Svings skarpt til høyre mot {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving",
                "name": "Ta en U-sving til {way_name}",
                "destination": "Ta en U-sving mot {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Hold {modifier} kjørefelt",
                "name": "Hold {modifier} kjørefelt inn på {way_name}",
                "destination": "Hold {modifier} kjørefelt mot {destination}"
            },
            "straight": {
                "default": "Hold kjørefelt",
                "name": "Hold kjørefelt inn på {way_name}",
                "destination": "Hold kjørefelt mot {destination}"
            },
            "slight left": {
                "default": "Hold venstre kjørefelt",
                "name": "Hold venstre kjørefelt inn på {way_name}",
                "destination": "Hold venstre kjørefelt mot {destination}"
            },
            "slight right": {
                "default": "Hold høyre kjørefelt",
                "name": "Hold høyre kjørefelt inn på {way_name}",
                "destination": "Hold høyre kjørefelt mot {destination}"
            },
            "sharp left": {
                "default": "Hold venstre kjørefelt",
                "name": "Hold venstre kjørefelt inn på {way_name}",
                "destination": "Hold venstre kjørefelt mot {destination}"
            },
            "sharp right": {
                "default": "Hold høyre kjørefelt",
                "name": "Hold høyre kjørefelt inn på {way_name}",
                "destination": "Hold høyre kjørefelt mot {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving",
                "name": "Ta en U-sving til {way_name}",
                "destination": "Ta en U-sving mot {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Fortsett {modifier}",
                "name": "Fortsett {modifier} til {way_name}",
                "destination": "Fortsett {modifier} mot  {destination}"
            },
            "straight": {
                "default": "Fortsett rett frem",
                "name": "Fortsett inn på {way_name}",
                "destination": "Fortsett mot {destination}"
            },
            "sharp left": {
                "default": "Sving skarpt til venstre",
                "name": "Sving skarpt til venstre inn på {way_name}",
                "destination": "Sving skarpt til venstre mot {destination}"
            },
            "sharp right": {
                "default": "Sving skarpt til høyre",
                "name": "Sving skarpt til høyre inn på {way_name}",
                "destination": "Svings skarpt til høyre mot {destination}"
            },
            "slight left": {
                "default": "Fortsett litt mot venstre",
                "name": "Fortsett litt mot venstre til {way_name}",
                "destination": "Fortsett litt mot venstre mot {destination}"
            },
            "slight right": {
                "default": "Fortsett litt mot høyre",
                "name": "Fortsett litt mot høyre til {way_name}",
                "destination": "Fortsett litt mot høyre mot {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving",
                "name": "Ta en U-sving til {way_name}",
                "destination": "Ta en U-sving mot {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Fortsett {modifier}",
                "name": "Fortsett {modifier} til {way_name}",
                "destination": "Fortsett {modifier} mot  {destination}"
            },
            "uturn": {
                "default": "Ta en U-sving",
                "name": "Ta en U-sving til {way_name}",
                "destination": "Ta en U-sving mot {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Ta avkjørselen",
                "name": "Ta avkjørselen inn på {way_name}",
                "destination": "Ta avkjørselen mot {destination}",
                "exit": "Ta avkjørsel {exit}",
                "exit_destination": "Ta avkjørsel {exit} mot {destination}"
            },
            "left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på venstre side",
                "exit_destination": "Ta avkjørsel {exit} på venstre side mot {destination}"
            },
            "right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på høyre side",
                "exit_destination": "Ta avkjørsel {exit} på høyre side mot {destination}"
            },
            "sharp left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på venstre side",
                "exit_destination": "Ta avkjørsel {exit} på venstre side mot {destination}"
            },
            "sharp right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på høyre side",
                "exit_destination": "Ta avkjørsel {exit} på høyre side mot {destination}"
            },
            "slight left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på venstre side",
                "exit_destination": "Ta avkjørsel {exit} på venstre side mot {destination}"
            },
            "slight right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}",
                "exit": "Ta avkjørsel {exit} på høyre side",
                "exit_destination": "Ta avkjørsel {exit} på høyre side mot {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Ta avkjørselen",
                "name": "Ta avkjørselen inn på {way_name}",
                "destination": "Ta avkjørselen mot {destination}"
            },
            "left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}"
            },
            "right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}"
            },
            "sharp left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}"
            },
            "sharp right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}"
            },
            "slight left": {
                "default": "Ta avkjørselen på venstre side",
                "name": "Ta avkjørselen på venstre side inn på {way_name}",
                "destination": "Ta avkjørselen på venstre side mot {destination}"
            },
            "slight right": {
                "default": "Ta avkjørselen på høyre side",
                "name": "Ta avkjørselen på høyre side inn på {way_name}",
                "destination": "Ta avkjørselen på høyre side mot {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Kjør inn i rundkjøringen",
                    "name": "Kjør inn i rundkjøringen og deretter ut på {way_name}",
                    "destination": "Kjør inn i rundkjøringen og deretter ut mot {destination}"
                },
                "name": {
                    "default": "Kjør inn i {rotary_name}",
                    "name": "Kjør inn i {rotary_name} og deretter ut på {way_name}",
                    "destination": "Kjør inn i {rotary_name} og deretter ut mot {destination}"
                },
                "exit": {
                    "default": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel",
                    "name": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel ut på {way_name}",
                    "destination": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel ut mot {destination} "
                },
                "name_exit": {
                    "default": "Kjør inn i {rotary_name} og ta {exit_number} avkjørsel",
                    "name": "Kjør inn i {rotary_name} og ta {exit_number} avkjørsel inn på {way_name}",
                    "destination": "Kjør inn i {rotary_name} og ta {exit_number} avkjørsel mot {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel",
                    "name": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel inn på {way_name}",
                    "destination": "Kjør inn i rundkjøringen og ta {exit_number} avkjørsel ut mot {destination} "
                },
                "default": {
                    "default": "Kjør inn i rundkjøringen",
                    "name": "Kjør inn i rundkjøringen og deretter ut på {way_name}",
                    "destination": "Kjør inn i rundkjøringen og deretter ut mot {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Ta en {modifier}",
                "name": "Ta en {modifier} inn på {way_name}",
                "destination": "Ta en {modifier} mot {destination}"
            },
            "left": {
                "default": "Sving til venstre",
                "name": "Sving til venstre inn på {way_name}",
                "destination": "Sving til venstre mot {destination}"
            },
            "right": {
                "default": "Sving til høyre",
                "name": "Sving til høyre inn på {way_name}",
                "destination": "Sving til høyre mot {destination}"
            },
            "straight": {
                "default": "Fortsett rett frem",
                "name": "Fortsett rett frem til  {way_name}",
                "destination": "Fortsett rett frem mot {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Kjør ut av rundkjøringen",
                "name": "Kjør ut av rundkjøringen og inn på {way_name}",
                "destination": "Kjør ut av rundkjøringen mot {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Kjør ut av rundkjøringen",
                "name": "Kjør ut av rundkjøringen og inn på {way_name}",
                "destination": "Kjør ut av rundkjøringen mot {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Ta en {modifier}",
                "name": "Ta en {modifier} inn på {way_name}",
                "destination": "Ta en {modifier} mot {destination}"
            },
            "left": {
                "default": "Sving til venstre",
                "name": "Sving til venstre inn på {way_name}",
                "destination": "Sving til venstre mot {destination}"
            },
            "right": {
                "default": "Sving til høyre",
                "name": "Sving til høyre inn på {way_name}",
                "destination": "Sving til høyre mot {destination}"
            },
            "straight": {
                "default": "Kjør rett frem",
                "name": "Kjør rett frem og inn på {way_name}",
                "destination": "Kjør rett frem mot {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Fortsett rett frem"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],38:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1.",
                "2": "2.",
                "3": "3.",
                "4": "4.",
                "5": "5.",
                "6": "6.",
                "7": "7.",
                "8": "8.",
                "9": "9.",
                "10": "10."
            },
            "direction": {
                "north": "północ",
                "northeast": "północny wschód",
                "east": "wschód",
                "southeast": "południowy wschód",
                "south": "południe",
                "southwest": "południowy zachód",
                "west": "zachód",
                "northwest": "północny zachód"
            },
            "modifier": {
                "left": "lewo",
                "right": "prawo",
                "sharp left": "ostro w lewo",
                "sharp right": "ostro w prawo",
                "slight left": "łagodnie w lewo",
                "slight right": "łagodnie w prawo",
                "straight": "prosto",
                "uturn": "zawróć"
            },
            "lanes": {
                "xo": "Trzymaj się prawej strony",
                "ox": "Trzymaj się lewej strony",
                "xox": "Trzymaj się środka",
                "oxo": "Trzymaj się lewej lub prawej strony"
            }
        },
        "modes": {
            "ferry": {
                "default": "Weź prom",
                "name": "Weź prom {way_name}",
                "destination": "Weź prom w kierunku {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, następnie za {distance} {instruction_two}",
            "two linked": "{instruction_one}, następnie {instruction_two}",
            "one in distance": "Za {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "Dojechano do miejsca docelowego {nth}",
                "upcoming": "Dojechano do miejsca docelowego {nth}",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}"
            },
            "left": {
                "default": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po lewej stronie"
            },
            "right": {
                "default": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po prawej stronie"
            },
            "sharp left": {
                "default": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po lewej stronie"
            },
            "sharp right": {
                "default": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po prawej stronie"
            },
            "slight right": {
                "default": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po prawej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po prawej stronie"
            },
            "slight left": {
                "default": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "upcoming": "Dojechano do miejsca docelowego {nth}, po lewej stronie",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, po lewej stronie"
            },
            "straight": {
                "default": "Dojechano do miejsca docelowego {nth} , prosto",
                "upcoming": "Dojechano do miejsca docelowego {nth} , prosto",
                "short": "Dojechano do miejsca docelowego {nth}",
                "short-upcoming": "Dojechano do miejsca docelowego {nth}",
                "named": "Dojechano do {waypoint_name}, prosto"
            }
        },
        "continue": {
            "default": {
                "default": "Skręć {modifier}",
                "name": "Skręć w {modifier}, aby pozostać na {way_name}",
                "destination": "Skręć {modifier} w kierunku {destination}",
                "exit": "Skręć {modifier} na {way_name}"
            },
            "straight": {
                "default": "Kontynuuj prosto",
                "name": "Jedź dalej prosto, aby pozostać na {way_name}",
                "destination": "Kontynuuj w kierunku {destination}",
                "distance": "Jedź dalej prosto przez {distance}",
                "namedistance": "Jedź dalej {way_name} przez {distance}"
            },
            "sharp left": {
                "default": "Skręć ostro w lewo",
                "name": "Skręć w lewo w ostry zakręt, aby pozostać na {way_name}",
                "destination": "Skręć ostro w lewo w kierunku {destination}"
            },
            "sharp right": {
                "default": "Skręć ostro w prawo",
                "name": "Skręć w prawo w ostry zakręt, aby pozostać na {way_name}",
                "destination": "Skręć ostro w prawo w kierunku {destination}"
            },
            "slight left": {
                "default": "Skręć w lewo w łagodny zakręt",
                "name": "Skręć w lewo w łagodny zakręt, aby pozostać na {way_name}",
                "destination": "Skręć w lewo w łagodny zakręt na {destination}"
            },
            "slight right": {
                "default": "Skręć w prawo w łagodny zakręt",
                "name": "Skręć w prawo w łagodny zakręt, aby pozostać na {way_name}",
                "destination": "Skręć w prawo w łagodny zakręt na {destination}"
            },
            "uturn": {
                "default": "Zawróć",
                "name": "Zawróć i jedź dalej {way_name}",
                "destination": "Zawróć w kierunku {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Kieruj się {direction}",
                "name": "Kieruj się {direction} na {way_name}",
                "namedistance": "Head {direction} on {way_name} for {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Skręć {modifier}",
                "name": "Skręć {modifier} na {way_name}",
                "destination": "Skręć {modifier} w kierunku {destination}"
            },
            "straight": {
                "default": "Kontynuuj prosto",
                "name": "Kontynuuj prosto na {way_name}",
                "destination": "Kontynuuj prosto w kierunku {destination}"
            },
            "uturn": {
                "default": "Zawróć na końcu ulicy",
                "name": "Zawróć na końcu ulicy na {way_name}",
                "destination": "Zawróć na końcu ulicy w kierunku {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "Na rozwidleniu trzymaj się {modifier}",
                "name": "Na rozwidleniu trzymaj się {modifier} na {way_name}",
                "destination": "Na rozwidleniu trzymaj się {modifier} w kierunku {destination}"
            },
            "slight left": {
                "default": "Na rozwidleniu trzymaj się lewej strony",
                "name": "Na rozwidleniu trzymaj się lewej strony w {way_name}",
                "destination": "Na rozwidleniu trzymaj się lewej strony w kierunku {destination}"
            },
            "slight right": {
                "default": "Na rozwidleniu trzymaj się prawej strony",
                "name": "Na rozwidleniu trzymaj się prawej strony na {way_name}",
                "destination": "Na rozwidleniu trzymaj się prawej strony w kierunku {destination}"
            },
            "sharp left": {
                "default": "Na rozwidleniu skręć ostro w lewo",
                "name": "Skręć ostro w lewo w {way_name}",
                "destination": "Skręć ostro w lewo w kierunku {destination}"
            },
            "sharp right": {
                "default": "Na rozwidleniu skręć ostro w prawo",
                "name": "Skręć ostro w prawo na {way_name}",
                "destination": "Skręć ostro w prawo w kierunku {destination}"
            },
            "uturn": {
                "default": "Zawróć",
                "name": "Zawróć na {way_name}",
                "destination": "Zawróć w kierunku {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Włącz się {modifier}",
                "name": "Włącz się {modifier} na {way_name}",
                "destination": "Włącz się {modifier} w kierunku {destination}"
            },
            "straight": {
                "default": "Włącz się prosto",
                "name": "Włącz się prosto na {way_name}",
                "destination": "Włącz się prosto w kierunku {destination}"
            },
            "slight left": {
                "default": "Włącz się z lewej strony",
                "name": "Włącz się z lewej strony na {way_name}",
                "destination": "Włącz się z lewej strony w kierunku {destination}"
            },
            "slight right": {
                "default": "Włącz się z prawej strony",
                "name": "Włącz się z prawej strony na {way_name}",
                "destination": "Włącz się z prawej strony w kierunku {destination}"
            },
            "sharp left": {
                "default": "Włącz się z lewej strony",
                "name": "Włącz się z lewej strony na {way_name}",
                "destination": "Włącz się z lewej strony w kierunku {destination}"
            },
            "sharp right": {
                "default": "Włącz się z prawej strony",
                "name": "Włącz się z prawej strony na {way_name}",
                "destination": "Włącz się z prawej strony w kierunku {destination}"
            },
            "uturn": {
                "default": "Zawróć",
                "name": "Zawróć na {way_name}",
                "destination": "Zawróć w kierunku {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Kontynuuj {modifier}",
                "name": "Kontynuuj {modifier} na {way_name}",
                "destination": "Kontynuuj {modifier} w kierunku {destination}"
            },
            "straight": {
                "default": "Kontynuuj prosto",
                "name": "Kontynuuj na {way_name}",
                "destination": "Kontynuuj w kierunku {destination}"
            },
            "sharp left": {
                "default": "Skręć ostro w lewo",
                "name": "Skręć ostro w lewo w {way_name}",
                "destination": "Skręć ostro w lewo w kierunku {destination}"
            },
            "sharp right": {
                "default": "Skręć ostro w prawo",
                "name": "Skręć ostro w prawo na {way_name}",
                "destination": "Skręć ostro w prawo w kierunku {destination}"
            },
            "slight left": {
                "default": "Kontynuuj łagodnie w lewo",
                "name": "Kontynuuj łagodnie w lewo na {way_name}",
                "destination": "Kontynuuj łagodnie w lewo w kierunku {destination}"
            },
            "slight right": {
                "default": "Kontynuuj łagodnie w prawo",
                "name": "Kontynuuj łagodnie w prawo na {way_name}",
                "destination": "Kontynuuj łagodnie w prawo w kierunku {destination}"
            },
            "uturn": {
                "default": "Zawróć",
                "name": "Zawróć na {way_name}",
                "destination": "Zawróć w kierunku {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Kontynuuj {modifier}",
                "name": "Kontynuuj {modifier} na {way_name}",
                "destination": "Kontynuuj {modifier} w kierunku {destination}"
            },
            "uturn": {
                "default": "Zawróć",
                "name": "Zawróć na {way_name}",
                "destination": "Zawróć w kierunku {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Zjedź",
                "name": "Weź zjazd na {way_name}",
                "destination": "Weź zjazd w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit}",
                "exit_destination": "Zjedź zjazdem {exit} na {destination}"
            },
            "left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po lewej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po lewej stronie na {destination}"
            },
            "right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po prawej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po prawej stronie na {destination}"
            },
            "sharp left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po lewej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po lewej stronie na {destination}"
            },
            "sharp right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po prawej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po prawej stronie na {destination}"
            },
            "slight left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po lewej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po lewej stronie na {destination}"
            },
            "slight right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}",
                "exit": "Zjedź zjazdem {exit} po prawej stronie",
                "exit_destination": "Zjedź zjazdem {exit} po prawej stronie na {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Weź zjazd",
                "name": "Weź zjazd na {way_name}",
                "destination": "Weź zjazd w kierunku {destination}"
            },
            "left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}"
            },
            "right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}"
            },
            "sharp left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}"
            },
            "sharp right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}"
            },
            "slight left": {
                "default": "Weź zjazd po lewej",
                "name": "Weź zjazd po lewej na {way_name}",
                "destination": "Weź zjazd po lewej w kierunku {destination}"
            },
            "slight right": {
                "default": "Weź zjazd po prawej",
                "name": "Weź zjazd po prawej na {way_name}",
                "destination": "Weź zjazd po prawej w kierunku {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Wjedź na rondo",
                    "name": "Wjedź na rondo i skręć na {way_name}",
                    "destination": "Wjedź na rondo i skręć w kierunku {destination}"
                },
                "name": {
                    "default": "Wjedź na {rotary_name}",
                    "name": "Wjedź na {rotary_name} i skręć na {way_name}",
                    "destination": "Wjedź na {rotary_name} i skręć w kierunku {destination}"
                },
                "exit": {
                    "default": "Wjedź na rondo i wyjedź {exit_number} zjazdem",
                    "name": "Wjedź na rondo i wyjedź {exit_number} zjazdem na {way_name}",
                    "destination": "Wjedź na rondo i wyjedź {exit_number} zjazdem w kierunku {destination}"
                },
                "name_exit": {
                    "default": "Wjedź na {rotary_name} i wyjedź {exit_number} zjazdem",
                    "name": "Wjedź na {rotary_name} i wyjedź {exit_number} zjazdem na {way_name}",
                    "destination": "Wjedź na {rotary_name} i wyjedź {exit_number} zjazdem w kierunku {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Wjedź na rondo i wyjedź {exit_number} zjazdem",
                    "name": "Wjedź na rondo i wyjedź {exit_number} zjazdem na {way_name}",
                    "destination": "Wjedź na rondo i wyjedź {exit_number} zjazdem w kierunku {destination}"
                },
                "default": {
                    "default": "Wjedź na rondo",
                    "name": "Wjedź na rondo i wyjedź na {way_name}",
                    "destination": "Wjedź na rondo i wyjedź w kierunku {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier}",
                "name": "{modifier} na {way_name}",
                "destination": "{modifier} w kierunku {destination}"
            },
            "left": {
                "default": "Skręć w lewo",
                "name": "Skręć w lewo na {way_name}",
                "destination": "Skręć w lewo w kierunku {destination}"
            },
            "right": {
                "default": "Skręć w prawo",
                "name": "Skręć w prawo na {way_name}",
                "destination": "Skręć w prawo w kierunku {destination}"
            },
            "straight": {
                "default": "Kontynuuj prosto",
                "name": "Kontynuuj prosto na {way_name}",
                "destination": "Kontynuuj prosto w kierunku {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "{modifier}",
                "name": "{modifier} na {way_name}",
                "destination": "{modifier} w kierunku {destination}"
            },
            "left": {
                "default": "Skręć w lewo",
                "name": "Skręć w lewo na {way_name}",
                "destination": "Skręć w lewo w kierunku {destination}"
            },
            "right": {
                "default": "Skręć w prawo",
                "name": "Skręć w prawo na {way_name}",
                "destination": "Skręć w prawo w kierunku {destination}"
            },
            "straight": {
                "default": "Kontynuuj prosto",
                "name": "Kontynuuj prosto na {way_name}",
                "destination": "Kontynuuj prosto w kierunku {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "{modifier}",
                "name": "{modifier} na {way_name}",
                "destination": "{modifier} w kierunku {destination}"
            },
            "left": {
                "default": "Skręć w lewo",
                "name": "Skręć w lewo na {way_name}",
                "destination": "Skręć w lewo w kierunku {destination}"
            },
            "right": {
                "default": "Skręć w prawo",
                "name": "Skręć w prawo na {way_name}",
                "destination": "Skręć w prawo w kierunku {destination}"
            },
            "straight": {
                "default": "Jedź prosto",
                "name": "Jedź prosto na {way_name}",
                "destination": "Jedź prosto w kierunku {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "{modifier}",
                "name": "{modifier} na {way_name}",
                "destination": "{modifier} w kierunku {destination}"
            },
            "left": {
                "default": "Skręć w lewo",
                "name": "Skręć w lewo na {way_name}",
                "destination": "Skręć w lewo w kierunku {destination}"
            },
            "right": {
                "default": "Skręć w prawo",
                "name": "Skręć w prawo na {way_name}",
                "destination": "Skręć w prawo w kierunku {destination}"
            },
            "straight": {
                "default": "Jedź prosto",
                "name": "Jedź prosto na {way_name}",
                "destination": "Jedź prosto w kierunku {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Kontynuuj prosto"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],39:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1º",
                "2": "2º",
                "3": "3º",
                "4": "4º",
                "5": "5º",
                "6": "6º",
                "7": "7º",
                "8": "8º",
                "9": "9º",
                "10": "10º"
            },
            "direction": {
                "north": "norte",
                "northeast": "nordeste",
                "east": "leste",
                "southeast": "sudeste",
                "south": "sul",
                "southwest": "sudoeste",
                "west": "oeste",
                "northwest": "noroeste"
            },
            "modifier": {
                "left": "à esquerda",
                "right": "à direita",
                "sharp left": "fechada à esquerda",
                "sharp right": "fechada à direita",
                "slight left": "suave à esquerda",
                "slight right": "suave à direita",
                "straight": "em frente",
                "uturn": "retorno"
            },
            "lanes": {
                "xo": "Mantenha-se à direita",
                "ox": "Mantenha-se à esquerda",
                "xox": "Mantenha-se ao centro",
                "oxo": "Mantenha-se à esquerda ou direita"
            }
        },
        "modes": {
            "ferry": {
                "default": "Pegue a balsa",
                "name": "Pegue a balsa {way_name}",
                "destination": "Pegue a balsa sentido {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, então, em {distance}, {instruction_two}",
            "two linked": "{instruction_one}, então {instruction_two}",
            "one in distance": "Em {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "saída {exit}"
        },
        "arrive": {
            "default": {
                "default": "Você chegou ao seu {nth} destino",
                "upcoming": "Você chegará ao seu {nth} destino",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou a {waypoint_name}"
            },
            "left": {
                "default": "Você chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Você chegará ao seu {nth} destino, à esquerda",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à esquerda"
            },
            "right": {
                "default": "Você chegou ao seu {nth} destino, à direita",
                "upcoming": "Você chegará ao seu {nth} destino, à direita",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à direita"
            },
            "sharp left": {
                "default": "Você chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Você chegará ao seu {nth} destino, à esquerda",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à esquerda"
            },
            "sharp right": {
                "default": "Você chegou ao seu {nth} destino, à direita",
                "upcoming": "Você chegará ao seu {nth} destino, à direita",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à direita"
            },
            "slight right": {
                "default": "Você chegou ao seu {nth} destino, à direita",
                "upcoming": "Você chegará ao seu {nth} destino, à direita",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à direita"
            },
            "slight left": {
                "default": "Você chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Você chegará ao seu {nth} destino, à esquerda",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "Você chegou {waypoint_name}, à esquerda"
            },
            "straight": {
                "default": "Você chegou ao seu {nth} destino, em frente",
                "upcoming": "Você vai chegar ao seu {nth} destino, em frente",
                "short": "Você chegou",
                "short-upcoming": "Você vai chegar",
                "named": "You have arrived at {waypoint_name}, straight ahead"
            }
        },
        "continue": {
            "default": {
                "default": "Vire {modifier}",
                "name": "Vire {modifier} para manter-se na {way_name}",
                "destination": "Vire {modifier} sentido {destination}",
                "exit": "Vire {modifier} em {way_name}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente para manter-se na {way_name}",
                "destination": "Continue em direção à {destination}",
                "distance": "Continue em frente por {distance}",
                "namedistance": "Continue na {way_name} por {distance}"
            },
            "sharp left": {
                "default": "Faça uma curva fechada a esquerda",
                "name": "Faça uma curva fechada a esquerda para manter-se na {way_name}",
                "destination": "Faça uma curva fechada a esquerda sentido {destination}"
            },
            "sharp right": {
                "default": "Faça uma curva fechada a direita",
                "name": "Faça uma curva fechada a direita para manter-se na {way_name}",
                "destination": "Faça uma curva fechada a direita sentido {destination}"
            },
            "slight left": {
                "default": "Faça uma curva suave a esquerda",
                "name": "Faça uma curva suave a esquerda para manter-se na {way_name}",
                "destination": "Faça uma curva suave a esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Faça uma curva suave a direita",
                "name": "Faça uma curva suave a direita para manter-se na {way_name}",
                "destination": "Faça uma curva suave a direita em direção a {destination}"
            },
            "uturn": {
                "default": "Faça o retorno",
                "name": "Faça o retorno e continue em {way_name}",
                "destination": "Faça o retorno sentido {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Siga {direction}",
                "name": "Siga {direction} em {way_name}",
                "namedistance": "Siga {direction} na {way_name} por {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Vire {modifier}",
                "name": "Vire {modifier} em {way_name}",
                "destination": "Vire {modifier} sentido {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente em {way_name}",
                "destination": "Continue em frente sentido {destination}"
            },
            "uturn": {
                "default": "Faça o retorno no fim da rua",
                "name": "Faça o retorno em {way_name} no fim da rua",
                "destination": "Faça o retorno sentido {destination} no fim da rua"
            }
        },
        "fork": {
            "default": {
                "default": "Mantenha-se {modifier} na bifurcação",
                "name": "Mantenha-se {modifier} na bifurcação em {way_name}",
                "destination": "Mantenha-se {modifier} na bifurcação sentido {destination}"
            },
            "slight left": {
                "default": "Mantenha-se à esquerda na bifurcação",
                "name": "Mantenha-se à esquerda na bifurcação em {way_name}",
                "destination": "Mantenha-se à esquerda na bifurcação sentido {destination}"
            },
            "slight right": {
                "default": "Mantenha-se à direita na bifurcação",
                "name": "Mantenha-se à direita na bifurcação em {way_name}",
                "destination": "Mantenha-se à direita na bifurcação sentido {destination}"
            },
            "sharp left": {
                "default": "Faça uma curva fechada à esquerda na bifurcação",
                "name": "Faça uma curva fechada à esquerda em {way_name}",
                "destination": "Faça uma curva fechada à esquerda sentido {destination}"
            },
            "sharp right": {
                "default": "Faça uma curva fechada à direita na bifurcação",
                "name": "Faça uma curva fechada à direita em {way_name}",
                "destination": "Faça uma curva fechada à direita sentido {destination}"
            },
            "uturn": {
                "default": "Faça o retorno",
                "name": "Faça o retorno em {way_name}",
                "destination": "Faça o retorno sentido {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Entre {modifier}",
                "name": "Entre {modifier} na {way_name}",
                "destination": "Entre {modifier} em direção à {destination}"
            },
            "straight": {
                "default": "Mesclar",
                "name": "Entre reto na {way_name}",
                "destination": "Entre reto em direção à {destination}"
            },
            "slight left": {
                "default": "Entre à esquerda",
                "name": "Entre à esquerda na {way_name}",
                "destination": "Entre à esquerda em direção à {destination}"
            },
            "slight right": {
                "default": "Entre à direita",
                "name": "Entre à direita na {way_name}",
                "destination": "Entre à direita em direção à {destination}"
            },
            "sharp left": {
                "default": "Entre à esquerda",
                "name": "Entre à esquerda na {way_name}",
                "destination": "Entre à esquerda em direção à {destination}"
            },
            "sharp right": {
                "default": "Entre à direita",
                "name": "Entre à direita na {way_name}",
                "destination": "Entre à direita em direção à {destination}"
            },
            "uturn": {
                "default": "Faça o retorno",
                "name": "Faça o retorno em {way_name}",
                "destination": "Faça o retorno sentido {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} em {way_name}",
                "destination": "Continue {modifier} sentido {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em {way_name}",
                "destination": "Continue em direção à {destination}"
            },
            "sharp left": {
                "default": "Faça uma curva fechada à esquerda",
                "name": "Faça uma curva fechada à esquerda em {way_name}",
                "destination": "Faça uma curva fechada à esquerda sentido {destination}"
            },
            "sharp right": {
                "default": "Faça uma curva fechada à direita",
                "name": "Faça uma curva fechada à direita em {way_name}",
                "destination": "Faça uma curva fechada à direita sentido {destination}"
            },
            "slight left": {
                "default": "Continue ligeiramente à esquerda",
                "name": "Continue ligeiramente à esquerda em {way_name}",
                "destination": "Continue ligeiramente à esquerda sentido {destination}"
            },
            "slight right": {
                "default": "Continue ligeiramente à direita",
                "name": "Continue ligeiramente à direita em {way_name}",
                "destination": "Continue ligeiramente à direita sentido {destination}"
            },
            "uturn": {
                "default": "Faça o retorno",
                "name": "Faça o retorno em {way_name}",
                "destination": "Faça o retorno sentido {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} em {way_name}",
                "destination": "Continue {modifier} sentido {destination}"
            },
            "uturn": {
                "default": "Faça o retorno",
                "name": "Faça o retorno em {way_name}",
                "destination": "Faça o retorno sentido {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Pegue a rampa",
                "name": "Pegue a rampa em {way_name}",
                "destination": "Pegue a rampa sentido {destination}",
                "exit": "Pegue a saída {exit}",
                "exit_destination": "Pegue a saída {exit} em direção à {destination}"
            },
            "left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}",
                "exit": "Pegue a saída {exit} à esquerda",
                "exit_destination": "Pegue a saída {exit}  à esquerda em direção à {destination}"
            },
            "right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentido {destination}",
                "exit": "Pegue a saída {exit} à direita",
                "exit_destination": "Pegue a saída {exit} à direita em direção à {destination}"
            },
            "sharp left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}",
                "exit": "Pegue a saída {exit} à esquerda",
                "exit_destination": "Pegue a saída {exit}  à esquerda em direção à {destination}"
            },
            "sharp right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentido {destination}",
                "exit": "Pegue a saída {exit} à direita",
                "exit_destination": "Pegue a saída {exit} à direita em direção à {destination}"
            },
            "slight left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}",
                "exit": "Pegue a saída {exit} à esquerda",
                "exit_destination": "Pegue a saída {exit}  à esquerda em direção à {destination}"
            },
            "slight right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentido {destination}",
                "exit": "Pegue a saída {exit} à direita",
                "exit_destination": "Pegue a saída {exit} à direita em direção à {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Pegue a rampa",
                "name": "Pegue a rampa em {way_name}",
                "destination": "Pegue a rampa sentido {destination}"
            },
            "left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}"
            },
            "right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentid {destination}"
            },
            "sharp left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}"
            },
            "sharp right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentido {destination}"
            },
            "slight left": {
                "default": "Pegue a rampa à esquerda",
                "name": "Pegue a rampa à esquerda em {way_name}",
                "destination": "Pegue a rampa à esquerda sentido {destination}"
            },
            "slight right": {
                "default": "Pegue a rampa à direita",
                "name": "Pegue a rampa à direita em {way_name}",
                "destination": "Pegue a rampa à direita sentido {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Entre na rotatória",
                    "name": "Entre na rotatória e saia na {way_name}",
                    "destination": "Entre na rotatória e saia sentido {destination}"
                },
                "name": {
                    "default": "Entre em {rotary_name}",
                    "name": "Entre em {rotary_name} e saia em {way_name}",
                    "destination": "Entre em {rotary_name} e saia sentido {destination}"
                },
                "exit": {
                    "default": "Entre na rotatória e pegue a {exit_number} saída",
                    "name": "Entre na rotatória e pegue a {exit_number} saída na {way_name}",
                    "destination": "Entre na rotatória e pegue a {exit_number} saída sentido {destination}"
                },
                "name_exit": {
                    "default": "Entre em {rotary_name} e saia na {exit_number} saída",
                    "name": "Entre em {rotary_name} e saia na {exit_number} saída em {way_name}",
                    "destination": "Entre em {rotary_name} e saia na {exit_number} saída sentido {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Entre na rotatória e pegue a {exit_number} saída",
                    "name": "Entre na rotatória e pegue a {exit_number} saída na {way_name}",
                    "destination": "Entre na rotatória e pegue a {exit_number} saída sentido {destination}"
                },
                "default": {
                    "default": "Entre na rotatória",
                    "name": "Entre na rotatória e saia na {way_name}",
                    "destination": "Entre na rotatória e saia sentido {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Siga {modifier}",
                "name": "Siga {modifier} em {way_name}",
                "destination": "Siga {modifier} sentido {destination}"
            },
            "left": {
                "default": "Vire à esquerda",
                "name": "Vire à esquerda em {way_name}",
                "destination": "Vire à esquerda sentido {destination}"
            },
            "right": {
                "default": "Vire à direita",
                "name": "Vire à direita em {way_name}",
                "destination": "Vire à direita sentido {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente em {way_name}",
                "destination": "Continue em frente sentido {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Saia da rotatória",
                "name": "Exit the traffic circle onto {way_name}",
                "destination": "Exit the traffic circle towards {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Saia da rotatória",
                "name": "Exit the traffic circle onto {way_name}",
                "destination": "Exit the traffic circle towards {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Siga {modifier}",
                "name": "Siga {modifier} em {way_name}",
                "destination": "Siga {modifier} sentido {destination}"
            },
            "left": {
                "default": "Vire à esquerda",
                "name": "Vire à esquerda em {way_name}",
                "destination": "Vire à esquerda sentido {destination}"
            },
            "right": {
                "default": "Vire à direita",
                "name": "Vire à direita em {way_name}",
                "destination": "Vire à direita sentido {destination}"
            },
            "straight": {
                "default": "Siga em frente",
                "name": "Siga em frente em {way_name}",
                "destination": "Siga em frente sentido {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continue em frente"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],40:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1º",
                "2": "2º",
                "3": "3º",
                "4": "4º",
                "5": "5º",
                "6": "6º",
                "7": "7º",
                "8": "8º",
                "9": "9º",
                "10": "10º"
            },
            "direction": {
                "north": "norte",
                "northeast": "nordeste",
                "east": "este",
                "southeast": "sudeste",
                "south": "sul",
                "southwest": "sudoeste",
                "west": "oeste",
                "northwest": "noroeste"
            },
            "modifier": {
                "left": "à esquerda",
                "right": "à direita",
                "sharp left": "acentuadamente à esquerda",
                "sharp right": "acentuadamente à direita",
                "slight left": "ligeiramente à esquerda",
                "slight right": "ligeiramente à direita",
                "straight": "em frente",
                "uturn": "inversão de marcha"
            },
            "lanes": {
                "xo": "Mantenha-se à direita",
                "ox": "Mantenha-se à esquerda",
                "xox": "Mantenha-se ao meio",
                "oxo": "Mantenha-se à esquerda ou à direita"
            }
        },
        "modes": {
            "ferry": {
                "default": "Apanhe o ferry",
                "name": "Apanhe o ferry {way_name}",
                "destination": "Apanhe o ferry para {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, depois, a {distance}, {instruction_two}",
            "two linked": "{instruction_one}, depois {instruction_two}",
            "one in distance": "A {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "saída {exit}"
        },
        "arrive": {
            "default": {
                "default": "Chegou ao seu {nth} destino",
                "upcoming": "Está a chegar ao seu {nth} destino",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}"
            },
            "left": {
                "default": "Chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Está a chegar ao seu {nth} destino, à esquerda",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à esquerda"
            },
            "right": {
                "default": "Chegou ao seu {nth} destino, à direita",
                "upcoming": "Está a chegar ao seu {nth} destino, à direita",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à direita"
            },
            "sharp left": {
                "default": "Chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Está a chegar ao seu {nth} destino, à esquerda",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à esquerda"
            },
            "sharp right": {
                "default": "Chegou ao seu {nth} destino, à direita",
                "upcoming": "Está a chegar ao seu {nth} destino, à direita",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à direita"
            },
            "slight right": {
                "default": "Chegou ao seu {nth} destino, à direita",
                "upcoming": "Está a chegar ao seu {nth} destino, à direita",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à direita"
            },
            "slight left": {
                "default": "Chegou ao seu {nth} destino, à esquerda",
                "upcoming": "Está a chegar ao seu {nth} destino, à esquerda",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, à esquerda"
            },
            "straight": {
                "default": "Chegou ao seu {nth} destino, em frente",
                "upcoming": "Está a chegar ao seu {nth} destino, em frente",
                "short": "Chegou",
                "short-upcoming": "Está a chegar",
                "named": "Chegou a {waypoint_name}, em frente"
            }
        },
        "continue": {
            "default": {
                "default": "Vire {modifier}",
                "name": "Vire {modifier} para se manter em {way_name}",
                "destination": "Vire {modifier} em direção a {destination}",
                "exit": "Vire {modifier} para {way_name}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente para se manter em {way_name}",
                "destination": "Continue em direção a {destination}",
                "distance": "Continue em frente por {distance}",
                "namedistance": "Continue em {way_name} por {distance}"
            },
            "sharp left": {
                "default": "Vire acentuadamente à esquerda",
                "name": "Vire acentuadamente à esquerda para se manter em {way_name}",
                "destination": "Vire acentuadamente à esquerda em direção a {destination}"
            },
            "sharp right": {
                "default": "Vire acentuadamente à direita",
                "name": "Vire acentuadamente à direita para se manter em {way_name}",
                "destination": "Vire acentuadamente à direita em direção a {destination}"
            },
            "slight left": {
                "default": "Vire ligeiramente à esquerda",
                "name": "Vire ligeiramente à esquerda para se manter em {way_name}",
                "destination": "Vire ligeiramente à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Vire ligeiramente à direita",
                "name": "Vire ligeiramente à direita para se manter em {way_name}",
                "destination": "Vire ligeiramente à direita em direção a {destination}"
            },
            "uturn": {
                "default": "Faça inversão de marcha",
                "name": "Faça inversão de marcha e continue em {way_name}",
                "destination": "Faça inversão de marcha em direção a {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Dirija-se para {direction}",
                "name": "Dirija-se para {direction} em {way_name}",
                "namedistance": "Dirija-se para {direction} em {way_name} por {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Vire {modifier}",
                "name": "Vire {modifier} para {way_name}",
                "destination": "Vire {modifier} em direção a {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente para {way_name}",
                "destination": "Continue em frente em direção a {destination}"
            },
            "uturn": {
                "default": "No final da estrada faça uma inversão de marcha",
                "name": "No final da estrada faça uma inversão de marcha para {way_name} ",
                "destination": "No final da estrada faça uma inversão de marcha em direção a {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "Na bifurcação mantenha-se {modifier}",
                "name": "Mantenha-se {modifier} para {way_name}",
                "destination": "Mantenha-se {modifier} em direção a {destination}"
            },
            "slight left": {
                "default": "Na bifurcação mantenha-se à esquerda",
                "name": "Mantenha-se à esquerda para {way_name}",
                "destination": "Mantenha-se à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Na bifurcação mantenha-se à direita",
                "name": "Mantenha-se à direita para {way_name}",
                "destination": "Mantenha-se à direita em direção a {destination}"
            },
            "sharp left": {
                "default": "Na bifurcação vire acentuadamente à esquerda",
                "name": "Vire acentuadamente à esquerda para {way_name}",
                "destination": "Vire acentuadamente à esquerda em direção a {destination}"
            },
            "sharp right": {
                "default": "Na bifurcação vire acentuadamente à direita",
                "name": "Vire acentuadamente à direita para {way_name}",
                "destination": "Vire acentuadamente à direita em direção a {destination}"
            },
            "uturn": {
                "default": "Faça inversão de marcha",
                "name": "Faça inversão de marcha para {way_name}",
                "destination": "Faça inversão de marcha em direção a {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Una-se ao tráfego {modifier}",
                "name": "Una-se ao tráfego {modifier} para {way_name}",
                "destination": "Una-se ao tráfego {modifier} em direção a {destination}"
            },
            "straight": {
                "default": "Una-se ao tráfego",
                "name": " Una-se ao tráfego para {way_name}",
                "destination": "Una-se ao tráfego em direção a {destination}"
            },
            "slight left": {
                "default": "Una-se ao tráfego à esquerda",
                "name": "Una-se ao tráfego à esquerda para {way_name}",
                "destination": "Una-se ao tráfego à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Una-se ao tráfego à direita",
                "name": "Una-se ao tráfego à direita para {way_name}",
                "destination": "Una-se ao tráfego à direita em direção a {destination}"
            },
            "sharp left": {
                "default": "Una-se ao tráfego à esquerda",
                "name": "Una-se ao tráfego à esquerda para {way_name}",
                "destination": "Una-se ao tráfego à esquerda em direção a {destination}"
            },
            "sharp right": {
                "default": "Una-se ao tráfego à direita",
                "name": "Una-se ao tráfego à direita para {way_name}",
                "destination": "Una-se ao tráfego à direita em direção a {destination}"
            },
            "uturn": {
                "default": "Faça inversão de marcha",
                "name": "Faça inversão de marcha para {way_name}",
                "destination": "Faça inversão de marcha em direção a {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} para {way_name}",
                "destination": "Continue {modifier} em direção a {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue para {way_name}",
                "destination": "Continue em direção a {destination}"
            },
            "sharp left": {
                "default": "Vire acentuadamente à esquerda",
                "name": "Vire acentuadamente à esquerda para {way_name}",
                "destination": "Vire acentuadamente à esquerda em direção a{destination}"
            },
            "sharp right": {
                "default": "Vire acentuadamente à direita",
                "name": "Vire acentuadamente à direita para {way_name}",
                "destination": "Vire acentuadamente à direita em direção a {destination}"
            },
            "slight left": {
                "default": "Continue ligeiramente à esquerda",
                "name": "Continue ligeiramente à esquerda para {way_name}",
                "destination": "Continue ligeiramente à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Continue ligeiramente à direita",
                "name": "Continue ligeiramente à direita para {way_name}",
                "destination": "Continue ligeiramente à direita em direção a {destination}"
            },
            "uturn": {
                "default": "Faça inversão de marcha",
                "name": "Faça inversão de marcha para {way_name}",
                "destination": "Faça inversão de marcha em direção a {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continue {modifier}",
                "name": "Continue {modifier} para {way_name}",
                "destination": "Continue {modifier} em direção a {destination}"
            },
            "uturn": {
                "default": "Faça inversão de marcha",
                "name": "Faça inversão de marcha para {way_name}",
                "destination": "Faça inversão de marcha em direção a {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Saia na saída",
                "name": "Saia na saída para {way_name}",
                "destination": "Saia na saída em direção a {destination}",
                "exit": "Saia na saída {exit}",
                "exit_destination": "Saia na saída {exit} em direção a {destination}"
            },
            "left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}",
                "exit": "Saia na saída {exit} à esquerda",
                "exit_destination": "Saia na saída {exit} à esquerda em direção a {destination}"
            },
            "right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}",
                "exit": "Saia na saída {exit} à direita",
                "exit_destination": "Saia na saída {exit} à direita em direção a {destination}"
            },
            "sharp left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}",
                "exit": "Saia na saída {exit} à esquerda",
                "exit_destination": "Saia na saída {exit} à esquerda em direção a {destination}"
            },
            "sharp right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}",
                "exit": "Saia na saída {exit} à direita",
                "exit_destination": "Saia na saída {exit} à direita em direção a {destination}"
            },
            "slight left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}",
                "exit": "Saia na saída {exit} à esquerda",
                "exit_destination": "Saia na saída {exit} à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}",
                "exit": "Saia na saída {exit} à direita",
                "exit_destination": "Saia na saída {exit} à direita em direção a {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Saia na saída",
                "name": "Saia na saída para {way_name}",
                "destination": "Saia na saída em direção a {destination}"
            },
            "left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}"
            },
            "right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}"
            },
            "sharp left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}"
            },
            "sharp right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}"
            },
            "slight left": {
                "default": "Saia na saída à esquerda",
                "name": "Saia na saída à esquerda para {way_name}",
                "destination": "Saia na saída à esquerda em direção a {destination}"
            },
            "slight right": {
                "default": "Saia na saída à direita",
                "name": "Saia na saída à direita para {way_name}",
                "destination": "Saia na saída à direita em direção a {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Entre na rotunda",
                    "name": "Entre na rotunda e saia para {way_name}",
                    "destination": "Entre na rotunda e saia em direção a {destination}"
                },
                "name": {
                    "default": "Entre em {rotary_name}",
                    "name": "Entre em {rotary_name} e saia para {way_name}",
                    "destination": "Entre em {rotary_name} e saia em direção a {destination}"
                },
                "exit": {
                    "default": "Entre na rotunda e saia na saída {exit_number}",
                    "name": "Entre na rotunda e saia na saída {exit_number} para {way_name}",
                    "destination": "Entre na rotunda e saia na saída {exit_number} em direção a {destination}"
                },
                "name_exit": {
                    "default": "Entre em {rotary_name} e saia na saída {exit_number}",
                    "name": "Entre em {rotary_name} e saia na saída {exit_number} para {way_name}",
                    "destination": "Entre em{rotary_name} e saia na saída {exit_number} em direção a {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Entre na rotunda e saia na saída {exit_number}",
                    "name": "Entre na rotunda e saia na saída {exit_number} para {way_name}",
                    "destination": "Entre na rotunda e saia na saída {exit_number} em direção a {destination}"
                },
                "default": {
                    "default": "Entre na rotunda",
                    "name": "Entre na rotunda e saia para {way_name}",
                    "destination": "Entre na rotunda e saia em direção a {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Siga {modifier}",
                "name": "Siga {modifier} para {way_name}",
                "destination": "Siga {modifier} em direção a {destination}"
            },
            "left": {
                "default": "Vire à esquerda",
                "name": "Vire à esquerda para {way_name}",
                "destination": "Vire à esquerda em direção a {destination}"
            },
            "right": {
                "default": "Vire à direita",
                "name": "Vire à direita para {way_name}",
                "destination": "Vire à direita em direção a {destination}"
            },
            "straight": {
                "default": "Continue em frente",
                "name": "Continue em frente para {way_name}",
                "destination": "Continue em frente em direção a {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Saia da rotunda",
                "name": "Saia da rotunda para {way_name}",
                "destination": "Saia da rotunda em direção a {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Saia da rotunda",
                "name": "Saia da rotunda para {way_name}",
                "destination": "Saia da rotunda em direção a {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Siga {modifier}",
                "name": "Siga {modifier} para{way_name}",
                "destination": "Siga {modifier} em direção a {destination}"
            },
            "left": {
                "default": "Vire à esquerda",
                "name": "Vire à esquerda para {way_name}",
                "destination": "Vire à esquerda em direção a {destination}"
            },
            "right": {
                "default": "Vire à direita",
                "name": "Vire à direita para {way_name}",
                "destination": "Vire à direita em direção a {destination}"
            },
            "straight": {
                "default": "Vá em frente",
                "name": "Vá em frente para {way_name}",
                "destination": "Vá em frente em direção a {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Continue em frente"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],41:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "prima",
                "2": "a doua",
                "3": "a treia",
                "4": "a patra",
                "5": "a cincea",
                "6": "a șasea",
                "7": "a șaptea",
                "8": "a opta",
                "9": "a noua",
                "10": "a zecea"
            },
            "direction": {
                "north": "nord",
                "northeast": "nord-est",
                "east": "est",
                "southeast": "sud-est",
                "south": "sud",
                "southwest": "sud-vest",
                "west": "vest",
                "northwest": "nord-vest"
            },
            "modifier": {
                "left": "stânga",
                "right": "dreapta",
                "sharp left": "puternic stânga",
                "sharp right": "puternic dreapta",
                "slight left": "ușor stânga",
                "slight right": "ușor dreapta",
                "straight": "înainte",
                "uturn": "întoarcere"
            },
            "lanes": {
                "xo": "Țineți stânga",
                "ox": "Țineți dreapta",
                "xox": "Țineți pe mijloc",
                "oxo": "Țineți pe laterale"
            }
        },
        "modes": {
            "ferry": {
                "default": "Luați feribotul",
                "name": "Luați feribotul {way_name}",
                "destination": "Luați feribotul spre {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, apoi în {distance}, {instruction_two}",
            "two linked": "{instruction_one} apoi {instruction_two}",
            "one in distance": "În {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "ieșirea {exit}"
        },
        "arrive": {
            "default": {
                "default": "Ați ajuns la {nth} destinație",
                "upcoming": "Ați ajuns la {nth} destinație",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}"
            },
            "left": {
                "default": "Ați ajuns la {nth} destinație, pe stânga",
                "upcoming": "Ați ajuns la {nth} destinație, pe stânga",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe stânga"
            },
            "right": {
                "default": "Ați ajuns la {nth} destinație, pe dreapta",
                "upcoming": "Ați ajuns la {nth} destinație, pe dreapta",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe dreapta"
            },
            "sharp left": {
                "default": "Ați ajuns la {nth} destinație, pe stânga",
                "upcoming": "Ați ajuns la {nth} destinație, pe stânga",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe stânga"
            },
            "sharp right": {
                "default": "Ați ajuns la {nth} destinație, pe dreapta",
                "upcoming": "Ați ajuns la {nth} destinație, pe dreapta",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe dreapta"
            },
            "slight right": {
                "default": "Ați ajuns la {nth} destinație, pe dreapta",
                "upcoming": "Ați ajuns la {nth} destinație, pe dreapta",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe dreapta"
            },
            "slight left": {
                "default": "Ați ajuns la {nth} destinație, pe stânga",
                "upcoming": "Ați ajuns la {nth} destinație, pe stânga",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, pe stânga"
            },
            "straight": {
                "default": "Ați ajuns la {nth} destinație, în față",
                "upcoming": "Ați ajuns la {nth} destinație, în față",
                "short": "Ați ajuns",
                "short-upcoming": "Veți ajunge",
                "named": "Ați ajuns {waypoint_name}, în față"
            }
        },
        "continue": {
            "default": {
                "default": "Virați {modifier}",
                "name": "Virați {modifier} pe {way_name}",
                "destination": "Virați {modifier} spre {destination}",
                "exit": "Virați {modifier} pe {way_name}"
            },
            "straight": {
                "default": "Mergeți înainte",
                "name": "Mergeți înainte pe {way_name}",
                "destination": "Continuați spre {destination}",
                "distance": "Mergeți înainte pentru {distance}",
                "namedistance": "Continuați pe {way_name} pentru {distance}"
            },
            "sharp left": {
                "default": "Virați puternic la stânga",
                "name": "Virați puternic la stânga pe {way_name}",
                "destination": "Virați puternic la stânga spre {destination}"
            },
            "sharp right": {
                "default": "Virați puternic la dreapta",
                "name": "Virați puternic la dreapta pe {way_name}",
                "destination": "Virați puternic la dreapta spre {destination}"
            },
            "slight left": {
                "default": "Virați ușor la stânga",
                "name": "Virați ușor la stânga pe {way_name}",
                "destination": "Virați ușor la stânga spre {destination}"
            },
            "slight right": {
                "default": "Virați ușor la dreapta",
                "name": "Virați ușor la dreapta pe {way_name}",
                "destination": "Virați ușor la dreapta spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă",
                "name": "Întoarceți-vă și continuați pe {way_name}",
                "destination": "Întoarceți-vă spre {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Mergeți spre {direction}",
                "name": "Mergeți spre {direction} pe {way_name}",
                "namedistance": "Mergeți spre {direction} pe {way_name} pentru {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Virați {modifier}",
                "name": "Virați {modifier} pe {way_name}",
                "destination": "Virați {modifier} spre {destination}"
            },
            "straight": {
                "default": "Continuați înainte",
                "name": "Continuați înainte pe {way_name}",
                "destination": "Continuați înainte spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă la sfârșitul drumului",
                "name": "Întoarceți-vă pe {way_name} la sfârșitul drumului",
                "destination": "Întoarceți-vă spre {destination} la sfârșitul drumului"
            }
        },
        "fork": {
            "default": {
                "default": "Țineți {modifier} la bifurcație",
                "name": "Țineți {modifier} la bifurcație pe {way_name}",
                "destination": "Țineți {modifier} la bifurcație spre {destination}"
            },
            "slight left": {
                "default": "Țineți pe stânga la bifurcație",
                "name": "Țineți pe stânga la bifurcație pe {way_name}",
                "destination": "Țineți pe stânga la bifurcație spre {destination}"
            },
            "slight right": {
                "default": "Țineți pe dreapta la bifurcație",
                "name": "Țineți pe dreapta la bifurcație pe {way_name}",
                "destination": "Țineți pe dreapta la bifurcație spre {destination}"
            },
            "sharp left": {
                "default": "Virați puternic stânga la bifurcație",
                "name": "Virați puternic stânga la bifurcație pe {way_name}",
                "destination": "Virați puternic stânga la bifurcație spre {destination}"
            },
            "sharp right": {
                "default": "Virați puternic dreapta la bifurcație",
                "name": "Virați puternic dreapta la bifurcație pe {way_name}",
                "destination": "Virați puternic dreapta la bifurcație spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă",
                "name": "Întoarceți-vă pe {way_name}",
                "destination": "Întoarceți-vă spre {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Intrați în {modifier}",
                "name": "Intrați în {modifier} pe {way_name}",
                "destination": "Intrați în {modifier} spre {destination}"
            },
            "straight": {
                "default": "Intrați",
                "name": "Intrați pe {way_name}",
                "destination": "Intrați spre {destination}"
            },
            "slight left": {
                "default": "Intrați în stânga",
                "name": "Intrați în stânga pe {way_name}",
                "destination": "Intrați în stânga spre {destination}"
            },
            "slight right": {
                "default": "Intrați în dreapta",
                "name": "Intrați în dreapta pe {way_name}",
                "destination": "Intrați în dreapta spre {destination}"
            },
            "sharp left": {
                "default": "Intrați în stânga",
                "name": "Intrați în stânga pe {way_name}",
                "destination": "Intrați în stânga spre {destination}"
            },
            "sharp right": {
                "default": "Intrați în dreapta",
                "name": "Intrați în dreapta pe {way_name}",
                "destination": "Intrați în dreapta spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă",
                "name": "Întoarceți-vă pe {way_name}",
                "destination": "Întoarceți-vă spre {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Continuați {modifier}",
                "name": "Continuați {modifier} pe {way_name}",
                "destination": "Continuați {modifier} spre {destination}"
            },
            "straight": {
                "default": "Continuați înainte",
                "name": "Continuați pe {way_name}",
                "destination": "Continuați spre {destination}"
            },
            "sharp left": {
                "default": "Virați puternic la stânga",
                "name": "Virați puternic la stânga pe {way_name}",
                "destination": "Virați puternic la stânga spre {destination}"
            },
            "sharp right": {
                "default": "Virați puternic la dreapta",
                "name": "Virați puternic la dreapta pe {way_name}",
                "destination": "Virați puternic la dreapta spre {destination}"
            },
            "slight left": {
                "default": "Continuați ușor la stânga",
                "name": "Continuați ușor la stânga pe {way_name}",
                "destination": "Continuați ușor la stânga spre {destination}"
            },
            "slight right": {
                "default": "Continuați ușor la dreapta",
                "name": "Continuați ușor la dreapta pe {way_name}",
                "destination": "Continuați ușor la dreapta spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă",
                "name": "Întoarceți-vă pe {way_name}",
                "destination": "Întoarceți-vă spre {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Continuați {modifier}",
                "name": "Continuați {modifier} pe {way_name}",
                "destination": "Continuați {modifier} spre {destination}"
            },
            "uturn": {
                "default": "Întoarceți-vă",
                "name": "Întoarceți-vă pe {way_name}",
                "destination": "Întoarceți-vă spre {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Urmați breteaua",
                "name": "Urmați breteaua pe {way_name}",
                "destination": "Urmați breteaua spre {destination}",
                "exit": "Urmați ieșirea {exit}",
                "exit_destination": "Urmați ieșirea {exit} spre {destination}"
            },
            "left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}",
                "exit": "Urmați ieșirea {exit} pe stânga",
                "exit_destination": "Urmați ieșirea {exit} pe stânga spre {destination}"
            },
            "right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}",
                "exit": "Urmați ieșirea {exit} pe dreapta",
                "exit_destination": "Urmați ieșirea {exit} pe dreapta spre {destination}"
            },
            "sharp left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}",
                "exit": "Urmați ieșirea {exit} pe stânga",
                "exit_destination": "Urmați ieșirea {exit} pe stânga spre {destination}"
            },
            "sharp right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}",
                "exit": "Urmați ieșirea {exit} pe dreapta",
                "exit_destination": "Urmați ieșirea {exit} pe dreapta spre {destination}"
            },
            "slight left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}",
                "exit": "Urmați ieșirea {exit} pe stânga",
                "exit_destination": "Urmați ieșirea {exit} pe stânga spre {destination}"
            },
            "slight right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}",
                "exit": "Urmați ieșirea {exit} pe dreapta",
                "exit_destination": "Urmați ieșirea {exit} pe dreapta spre {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Urmați breteaua de intrare",
                "name": "Urmați breteaua pe {way_name}",
                "destination": "Urmați breteaua spre {destination}"
            },
            "left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}"
            },
            "right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}"
            },
            "sharp left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}"
            },
            "sharp right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}"
            },
            "slight left": {
                "default": "Urmați breteaua din stânga",
                "name": "Urmați breteaua din stânga pe {way_name}",
                "destination": "Urmați breteaua din stânga spre {destination}"
            },
            "slight right": {
                "default": "Urmați breteaua din dreapta",
                "name": "Urmați breteaua din dreapta pe {way_name}",
                "destination": "Urmați breteaua din dreapta spre {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Intrați în sensul giratoriu",
                    "name": "Intrați în sensul giratoriu și ieșiți pe {way_name}",
                    "destination": "Intrați în sensul giratoriu și ieșiți spre {destination}"
                },
                "name": {
                    "default": "Intrați în {rotary_name}",
                    "name": "Intrați în {rotary_name} și ieșiți pe {way_name}",
                    "destination": "Intrați în {rotary_name} și ieșiți spre {destination}"
                },
                "exit": {
                    "default": "Intrați în sensul giratoriu și urmați {exit_number} ieșire",
                    "name": "Intrați în sensul giratoriu și urmați {exit_number} ieșire pe {way_name}",
                    "destination": "Intrați în sensul giratoriu și urmați {exit_number} ieșire spre {destination}"
                },
                "name_exit": {
                    "default": "Intrați în {rotary_name} și urmați {exit_number} ieșire",
                    "name": "Intrați în {rotary_name} și urmați {exit_number} ieșire pe {way_name}",
                    "destination": "Intrați în  {rotary_name} și urmați {exit_number} ieșire spre {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Intrați în sensul giratoriu și urmați {exit_number} ieșire",
                    "name": "Intrați în sensul giratoriu și urmați {exit_number} ieșire pe {way_name}",
                    "destination": "Intrați în sensul giratoriu și urmați {exit_number} ieșire spre {destination}"
                },
                "default": {
                    "default": "Intrați în sensul giratoriu",
                    "name": "Intrați în sensul giratoriu și ieșiți pe {way_name}",
                    "destination": "Intrați în sensul giratoriu și ieșiți spre {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "La sensul giratoriu virați {modifier}",
                "name": "La sensul giratoriu virați {modifier} pe {way_name}",
                "destination": "La sensul giratoriu virați {modifier} spre {destination}"
            },
            "left": {
                "default": "La sensul giratoriu virați la stânga",
                "name": "La sensul giratoriu virați la stânga pe {way_name}",
                "destination": "La sensul giratoriu virați la stânga spre {destination}"
            },
            "right": {
                "default": "La sensul giratoriu virați la dreapta",
                "name": "La sensul giratoriu virați la dreapta pe {way_name}",
                "destination": "La sensul giratoriu virați la dreapta spre {destination}"
            },
            "straight": {
                "default": "La sensul giratoriu continuați înainte",
                "name": "La sensul giratoriu continuați înainte pe {way_name}",
                "destination": "La sensul giratoriu continuați înainte spre {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Ieșiți din sensul giratoriu",
                "name": "Ieșiți din sensul giratoriu pe {way_name}",
                "destination": "Ieșiți din sensul giratoriu spre {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Ieșiți din sensul giratoriu",
                "name": "Ieșiți din sensul giratoriu pe {way_name}",
                "destination": "Ieșiți din sensul giratoriu spre {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Virați {modifier}",
                "name": "Virați {modifier} pe {way_name}",
                "destination": "Virați {modifier} spre {destination}"
            },
            "left": {
                "default": "Virați la stânga",
                "name": "Virați la stânga pe {way_name}",
                "destination": "Virați la stânga spre {destination}"
            },
            "right": {
                "default": "Virați la dreapta",
                "name": "Virați la dreapta pe {way_name}",
                "destination": "Virați la dreapta spre {destination}"
            },
            "straight": {
                "default": "Mergeți înainte",
                "name": "Mergeți înainte pe {way_name}",
                "destination": "Mergeți înainte spre {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Mergeți înainte"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],42:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "первый",
                "2": "второй",
                "3": "третий",
                "4": "четвёртый",
                "5": "пятый",
                "6": "шестой",
                "7": "седьмой",
                "8": "восьмой",
                "9": "девятый",
                "10": "десятый"
            },
            "direction": {
                "north": "северном",
                "northeast": "северо-восточном",
                "east": "восточном",
                "southeast": "юго-восточном",
                "south": "южном",
                "southwest": "юго-западном",
                "west": "западном",
                "northwest": "северо-западном"
            },
            "modifier": {
                "left": "налево",
                "right": "направо",
                "sharp left": "налево",
                "sharp right": "направо",
                "slight left": "левее",
                "slight right": "правее",
                "straight": "прямо",
                "uturn": "на разворот"
            },
            "lanes": {
                "xo": "Держитесь правее",
                "ox": "Держитесь левее",
                "xox": "Держитесь посередине",
                "oxo": "Держитесь слева или справа"
            }
        },
        "modes": {
            "ferry": {
                "default": "Погрузитесь на паром",
                "name": "Погрузитесь на паром {way_name}",
                "destination": "Погрузитесь на паром в направлении {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, затем через {distance} {instruction_two}",
            "two linked": "{instruction_one}, затем {instruction_two}",
            "one in distance": "Через {distance} {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "съезд {exit}"
        },
        "arrive": {
            "default": {
                "default": "Вы прибыли в {nth} пункт назначения",
                "upcoming": "Вы прибудете в {nth} пункт назначения",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}"
            },
            "left": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится слева",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет слева",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится слева"
            },
            "right": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится справа",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет справа",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится справа"
            },
            "sharp left": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится слева сзади",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет слева сзади",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится слева сзади"
            },
            "sharp right": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится справа сзади",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет справа сзади",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится справа сзади"
            },
            "slight right": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится справа впереди",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет справа впереди",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится справа впереди"
            },
            "slight left": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится слева впереди",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет слева впереди",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится слева впереди"
            },
            "straight": {
                "default": "Вы прибыли в {nth} пункт назначения, он находится перед Вами",
                "upcoming": "Вы прибудете в {nth} пункт назначения, он будет перед Вами",
                "short": "Вы прибыли",
                "short-upcoming": "Вы скоро прибудете",
                "named": "Вы прибыли в пункт назначения, {waypoint_name}, он находится перед Вами"
            }
        },
        "continue": {
            "default": {
                "default": "Двигайтесь {modifier}",
                "name": "Двигайтесь {modifier} по {way_name:dative}",
                "destination": "Двигайтесь {modifier} в направлении {destination}",
                "exit": "Двигайтесь {modifier} на {way_name:accusative}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Продолжите движение по {way_name:dative}",
                "destination": "Продолжите движение в направлении {destination}",
                "distance": "Двигайтесь прямо {distance}",
                "namedistance": "Двигайтесь прямо {distance} по {way_name:dative}"
            },
            "sharp left": {
                "default": "Резко поверните налево",
                "name": "Резко поверните налево на {way_name:accusative}",
                "destination": "Резко поверните налево в направлении {destination}"
            },
            "sharp right": {
                "default": "Резко поверните направо",
                "name": "Резко поверните направо на {way_name:accusative}",
                "destination": "Резко поверните направо в направлении {destination}"
            },
            "slight left": {
                "default": "Плавно поверните налево",
                "name": "Плавно поверните налево на {way_name:accusative}",
                "destination": "Плавно поверните налево в направлении {destination}"
            },
            "slight right": {
                "default": "Плавно поверните направо",
                "name": "Плавно поверните направо на {way_name:accusative}",
                "destination": "Плавно поверните направо в направлении {destination}"
            },
            "uturn": {
                "default": "Развернитесь",
                "name": "Развернитесь и продолжите движение по {way_name:dative}",
                "destination": "Развернитесь в направлении {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Двигайтесь в {direction} направлении",
                "name": "Двигайтесь в {direction} направлении по {way_name:dative}",
                "namedistance": "Двигайтесь {distance} в {direction} направлении по {way_name:dative}"
            }
        },
        "end of road": {
            "default": {
                "default": "Поверните {modifier}",
                "name": "Поверните {modifier} на {way_name:accusative}",
                "destination": "Поверните {modifier} в направлении {destination}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Двигайтесь прямо по {way_name:dative}",
                "destination": "Двигайтесь прямо в направлении {destination}"
            },
            "uturn": {
                "default": "В конце дороги развернитесь",
                "name": "Развернитесь в конце {way_name:genitive}",
                "destination": "В конце дороги развернитесь в направлении {destination}"
            }
        },
        "fork": {
            "default": {
                "default": "На развилке двигайтесь {modifier}",
                "name": "На развилке двигайтесь {modifier} на {way_name:accusative}",
                "destination": "На развилке двигайтесь {modifier} в направлении {destination}"
            },
            "slight left": {
                "default": "На развилке держитесь левее",
                "name": "На развилке держитесь левее на {way_name:accusative}",
                "destination": "На развилке держитесь левее и продолжите движение в направлении {destination}"
            },
            "slight right": {
                "default": "На развилке держитесь правее",
                "name": "На развилке держитесь правее на {way_name:accusative}",
                "destination": "На развилке держитесь правее и продолжите движение в направлении {destination}"
            },
            "sharp left": {
                "default": "На развилке резко поверните налево",
                "name": "Резко поверните налево на {way_name:accusative}",
                "destination": "Резко поверните налево и продолжите движение в направлении {destination}"
            },
            "sharp right": {
                "default": "На развилке резко поверните направо",
                "name": "Резко поверните направо на {way_name:accusative}",
                "destination": "Резко поверните направо и продолжите движение в направлении {destination}"
            },
            "uturn": {
                "default": "На развилке развернитесь",
                "name": "На развилке развернитесь на {way_name:prepositional}",
                "destination": "На развилке развернитесь и продолжите движение в направлении {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Перестройтесь {modifier}",
                "name": "Перестройтесь {modifier} на {way_name:accusative}",
                "destination": "Перестройтесь {modifier} в направлении {destination}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Продолжите движение по {way_name:dative}",
                "destination": "Продолжите движение в направлении {destination}"
            },
            "slight left": {
                "default": "Перестройтесь левее",
                "name": "Перестройтесь левее на {way_name:accusative}",
                "destination": "Перестройтесь левее в направлении {destination}"
            },
            "slight right": {
                "default": "Перестройтесь правее",
                "name": "Перестройтесь правее на {way_name:accusative}",
                "destination": "Перестройтесь правее в направлении {destination}"
            },
            "sharp left": {
                "default": "Перестраивайтесь левее",
                "name": "Перестраивайтесь левее на {way_name:accusative}",
                "destination": "Перестраивайтесь левее в направлении {destination}"
            },
            "sharp right": {
                "default": "Перестраивайтесь правее",
                "name": "Перестраивайтесь правее на {way_name:accusative}",
                "destination": "Перестраивайтесь правее в направлении {destination}"
            },
            "uturn": {
                "default": "Развернитесь",
                "name": "Развернитесь на {way_name:prepositional}",
                "destination": "Развернитесь в направлении {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Двигайтесь {modifier}",
                "name": "Двигайтесь {modifier} на {way_name:accusative}",
                "destination": "Двигайтесь {modifier} в направлении {destination}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Продолжите движение по {way_name:dative}",
                "destination": "Продолжите движение в направлении {destination}"
            },
            "sharp left": {
                "default": "Резко поверните налево",
                "name": "Резко поверните налево на {way_name:accusative}",
                "destination": "Резко поверните налево и продолжите движение в направлении {destination}"
            },
            "sharp right": {
                "default": "Резко поверните направо",
                "name": "Резко поверните направо на {way_name:accusative}",
                "destination": "Резко поверните направо и продолжите движение в направлении {destination}"
            },
            "slight left": {
                "default": "Плавно поверните налево",
                "name": "Плавно поверните налево на {way_name:accusative}",
                "destination": "Плавно поверните налево в направлении {destination}"
            },
            "slight right": {
                "default": "Плавно поверните направо",
                "name": "Плавно поверните направо на {way_name:accusative}",
                "destination": "Плавно поверните направо в направлении {destination}"
            },
            "uturn": {
                "default": "Развернитесь",
                "name": "Развернитесь на {way_name:prepositional}",
                "destination": "Развернитесь и продолжите движение в направлении {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Двигайтесь {modifier}",
                "name": "Двигайтесь {modifier} по {way_name:dative}",
                "destination": "Двигайтесь {modifier} в направлении {destination}"
            },
            "uturn": {
                "default": "Развернитесь",
                "name": "Развернитесь на {way_name:prepositional}",
                "destination": "Развернитесь и продолжите движение в направлении {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Сверните на съезд",
                "name": "Сверните на съезд на {way_name:accusative}",
                "destination": "Сверните на съезд в направлении {destination}",
                "exit": "Сверните на съезд {exit}",
                "exit_destination": "Сверните на съезд {exit} в направлении {destination}"
            },
            "left": {
                "default": "Сверните на левый съезд",
                "name": "Сверните на левый съезд на {way_name:accusative}",
                "destination": "Сверните на левый съезд в направлении {destination}",
                "exit": "Сверните на съезд {exit} слева",
                "exit_destination": "Сверните на съезд {exit} слева в направлении {destination}"
            },
            "right": {
                "default": "Сверните на правый съезд",
                "name": "Сверните на правый съезд на {way_name:accusative}",
                "destination": "Сверните на правый съезд в направлении {destination}",
                "exit": "Сверните на съезд {exit} справа",
                "exit_destination": "Сверните на съезд {exit} справа в направлении {destination}"
            },
            "sharp left": {
                "default": "Поверните налево на съезд",
                "name": "Поверните налево на съезд на {way_name:accusative}",
                "destination": "Поверните налево на съезд в направлении {destination}",
                "exit": "Поверните налево на съезд {exit}",
                "exit_destination": "Поверните налево на съезд {exit} в направлении {destination}"
            },
            "sharp right": {
                "default": "Поверните направо на съезд",
                "name": "Поверните направо на съезд на {way_name:accusative}",
                "destination": "Поверните направо на съезд в направлении {destination}",
                "exit": "Поверните направо на съезд {exit}",
                "exit_destination": "Поверните направо на съезд {exit} в направлении {destination}"
            },
            "slight left": {
                "default": "Перестройтесь левее на съезд",
                "name": "Перестройтесь левее на съезд на {way_name:accusative}",
                "destination": "Перестройтесь левее на съезд в направлении {destination}",
                "exit": "Перестройтесь левее на {exit}",
                "exit_destination": "Перестройтесь левее на съезд {exit} в направлении {destination}"
            },
            "slight right": {
                "default": "Перестройтесь правее на съезд",
                "name": "Перестройтесь правее на съезд на {way_name:accusative}",
                "destination": "Перестройтесь правее на съезд в направлении {destination}",
                "exit": "Перестройтесь правее на съезд {exit}",
                "exit_destination": "Перестройтесь правее на съезд {exit} в направлении {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Сверните на автомагистраль",
                "name": "Сверните на въезд на {way_name:accusative}",
                "destination": "Сверните на въезд на автомагистраль в направлении {destination}"
            },
            "left": {
                "default": "Сверните на левый въезд на автомагистраль",
                "name": "Сверните на левый въезд на {way_name:accusative}",
                "destination": "Сверните на левый въезд на автомагистраль в направлении {destination}"
            },
            "right": {
                "default": "Сверните на правый въезд на автомагистраль",
                "name": "Сверните на правый въезд на {way_name:accusative}",
                "destination": "Сверните на правый въезд на автомагистраль в направлении {destination}"
            },
            "sharp left": {
                "default": "Поверните на левый въезд на автомагистраль",
                "name": "Поверните на левый въезд на {way_name:accusative}",
                "destination": "Поверните на левый въезд на автомагистраль в направлении {destination}"
            },
            "sharp right": {
                "default": "Поверните на правый въезд на автомагистраль",
                "name": "Поверните на правый въезд на {way_name:accusative}",
                "destination": "Поверните на правый въезд на автомагистраль в направлении {destination}"
            },
            "slight left": {
                "default": "Перестройтесь левее на въезд на автомагистраль",
                "name": "Перестройтесь левее на {way_name:accusative}",
                "destination": "Перестройтесь левее на автомагистраль в направлении {destination}"
            },
            "slight right": {
                "default": "Перестройтесь правее на въезд на автомагистраль",
                "name": "Перестройтесь правее на {way_name:accusative}",
                "destination": "Перестройтесь правее на автомагистраль в направлении {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Продолжите движение по круговой развязке",
                    "name": "На круговой развязке сверните на {way_name:accusative}",
                    "destination": "На круговой развязке сверните в направлении {destination}"
                },
                "name": {
                    "default": "Продолжите движение по {rotary_name:dative}",
                    "name": "На {rotary_name:prepositional} сверните на {way_name:accusative}",
                    "destination": "На {rotary_name:prepositional} сверните в направлении {destination}"
                },
                "exit": {
                    "default": "На круговой развязке сверните на {exit_number} съезд",
                    "name": "На круговой развязке сверните на {exit_number} съезд на {way_name:accusative}",
                    "destination": "На круговой развязке сверните на {exit_number} съезд в направлении {destination}"
                },
                "name_exit": {
                    "default": "На {rotary_name:prepositional} сверните на {exit_number} съезд",
                    "name": "На {rotary_name:prepositional} сверните на {exit_number} съезд на {way_name:accusative}",
                    "destination": "На {rotary_name:prepositional} сверните на {exit_number} съезд в направлении {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "На круговой развязке сверните на {exit_number} съезд",
                    "name": "На круговой развязке сверните на {exit_number} съезд на {way_name:accusative}",
                    "destination": "На круговой развязке сверните на {exit_number} съезд в направлении {destination}"
                },
                "default": {
                    "default": "Продолжите движение по круговой развязке",
                    "name": "На круговой развязке сверните на {way_name:accusative}",
                    "destination": "На круговой развязке сверните в направлении {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Двигайтесь {modifier}",
                "name": "Двигайтесь {modifier} на {way_name:accusative}",
                "destination": "Двигайтесь {modifier} в направлении {destination}"
            },
            "left": {
                "default": "Сверните налево",
                "name": "Сверните налево на {way_name:accusative}",
                "destination": "Сверните налево в направлении {destination}"
            },
            "right": {
                "default": "Сверните направо",
                "name": "Сверните направо на {way_name:accusative}",
                "destination": "Сверните направо в направлении {destination}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Двигайтесь прямо по {way_name:dative}",
                "destination": "Двигайтесь прямо в направлении {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Сверните с круговой развязки",
                "name": "Сверните с круговой развязки на {way_name:accusative}",
                "destination": "Сверните с круговой развязки в направлении {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Сверните с круговой развязки",
                "name": "Сверните с круговой развязки на {way_name:accusative}",
                "destination": "Сверните с круговой развязки в направлении {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Двигайтесь {modifier}",
                "name": "Двигайтесь {modifier} на {way_name:accusative}",
                "destination": "Двигайтесь {modifier}  в направлении {destination}"
            },
            "left": {
                "default": "Поверните налево",
                "name": "Поверните налево на {way_name:accusative}",
                "destination": "Поверните налево в направлении {destination}"
            },
            "right": {
                "default": "Поверните направо",
                "name": "Поверните направо на {way_name:accusative}",
                "destination": "Поверните направо  в направлении {destination}"
            },
            "straight": {
                "default": "Двигайтесь прямо",
                "name": "Двигайтесь по {way_name:dative}",
                "destination": "Двигайтесь в направлении {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Продолжайте движение прямо"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],43:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1:a",
                "2": "2:a",
                "3": "3:e",
                "4": "4:e",
                "5": "5:e",
                "6": "6:e",
                "7": "7:e",
                "8": "8:e",
                "9": "9:e",
                "10": "10:e"
            },
            "direction": {
                "north": "norr",
                "northeast": "nordost",
                "east": "öster",
                "southeast": "sydost",
                "south": "söder",
                "southwest": "sydväst",
                "west": "väster",
                "northwest": "nordväst"
            },
            "modifier": {
                "left": "vänster",
                "right": "höger",
                "sharp left": "vänster",
                "sharp right": "höger",
                "slight left": "vänster",
                "slight right": "höger",
                "straight": "rakt fram",
                "uturn": "U-sväng"
            },
            "lanes": {
                "xo": "Håll till höger",
                "ox": "Håll till vänster",
                "xox": "Håll till mitten",
                "oxo": "Håll till vänster eller höger"
            }
        },
        "modes": {
            "ferry": {
                "default": "Ta färjan",
                "name": "Ta färjan på {way_name}",
                "destination": "Ta färjan mot {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, sedan efter {distance}, {instruction_two}",
            "two linked": "{instruction_one}, sedan {instruction_two}",
            "one in distance": "Om {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "Du är framme vid din {nth} destination",
                "upcoming": "Du är snart framme vid din {nth} destination",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}"
            },
            "left": {
                "default": "Du är framme vid din {nth} destination, till vänster",
                "upcoming": "Du är snart framme vid din {nth} destination, till vänster",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till vänster"
            },
            "right": {
                "default": "Du är framme vid din {nth} destination, till höger",
                "upcoming": "Du är snart framme vid din {nth} destination, till höger",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till höger"
            },
            "sharp left": {
                "default": "Du är framme vid din {nth} destination, till vänster",
                "upcoming": "Du är snart framme vid din {nth} destination, till vänster",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till vänster"
            },
            "sharp right": {
                "default": "Du är framme vid din {nth} destination, till höger",
                "upcoming": "Du är snart framme vid din {nth} destination, till höger",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till höger"
            },
            "slight right": {
                "default": "Du är framme vid din {nth} destination, till höger",
                "upcoming": "Du är snart framme vid din {nth} destination, till höger",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till höger"
            },
            "slight left": {
                "default": "Du är framme vid din {nth} destination, till vänster",
                "upcoming": "Du är snart framme vid din {nth} destination, till vänster",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, till vänster"
            },
            "straight": {
                "default": "Du är framme vid din {nth} destination, rakt fram",
                "upcoming": "Du är snart framme vid din {nth} destination, rakt fram",
                "short": "Du är framme",
                "short-upcoming": "Du är snart framme",
                "named": "Du är framme vid {waypoint_name}, rakt fram"
            }
        },
        "continue": {
            "default": {
                "default": "Sväng {modifier}",
                "name": "Sväng {modifier} och fortsätt på {way_name}",
                "destination": "Sväng {modifier} mot {destination}",
                "exit": "Sväng {modifier} in på {way_name}"
            },
            "straight": {
                "default": "Fortsätt rakt fram",
                "name": "Kör rakt fram och fortsätt på {way_name}",
                "destination": "Fortsätt mot {destination}",
                "distance": "Fortsätt rakt fram i {distance}",
                "namedistance": "Fortsätt på {way_name} i {distance}"
            },
            "sharp left": {
                "default": "Sväng vänster",
                "name": "Sväng vänster och fortsätt på {way_name}",
                "destination": "Sväng vänster mot {destination}"
            },
            "sharp right": {
                "default": "Sväng höger",
                "name": "Sväng höger och fortsätt på {way_name}",
                "destination": "Sväng höger mot {destination}"
            },
            "slight left": {
                "default": "Sväng vänster",
                "name": "Sväng vänster och fortsätt på {way_name}",
                "destination": "Sväng vänster mot {destination}"
            },
            "slight right": {
                "default": "Sväng höger",
                "name": "Sväng höger och fortsätt på {way_name}",
                "destination": "Sväng höger mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng",
                "name": "Gör en U-sväng och fortsätt på {way_name}",
                "destination": "Gör en U-sväng mot {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Kör åt {direction}",
                "name": "Kör åt {direction} på {way_name}",
                "namedistance": "Kör {distance} åt {direction} på {way_name}"
            }
        },
        "end of road": {
            "default": {
                "default": "Sväng {modifier}",
                "name": "Sväng {modifier} in på {way_name}",
                "destination": "Sväng {modifier} mot {destination}"
            },
            "straight": {
                "default": "Fortsätt rakt fram",
                "name": "Fortsätt rakt fram in på {way_name}",
                "destination": "Fortsätt rakt fram mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng i slutet av vägen",
                "name": "Gör en U-sväng in på {way_name} i slutet av vägen",
                "destination": "Gör en U-sväng mot {destination} i slutet av vägen"
            }
        },
        "fork": {
            "default": {
                "default": "Håll till {modifier} där vägen delar sig",
                "name": "Håll till {modifier} in på {way_name}",
                "destination": "Håll till {modifier} mot {destination}"
            },
            "slight left": {
                "default": "Håll till vänster där vägen delar sig",
                "name": "Håll till vänster in på {way_name}",
                "destination": "Håll till vänster mot {destination}"
            },
            "slight right": {
                "default": "Håll till höger där vägen delar sig",
                "name": "Håll till höger in på {way_name}",
                "destination": "Håll till höger mot {destination}"
            },
            "sharp left": {
                "default": "Sväng vänster där vägen delar sig",
                "name": "Sväng vänster in på {way_name}",
                "destination": "Sväng vänster mot {destination}"
            },
            "sharp right": {
                "default": "Sväng höger där vägen delar sig",
                "name": "Sväng höger in på {way_name}",
                "destination": "Sväng höger mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng",
                "name": "Gör en U-sväng in på {way_name}",
                "destination": "Gör en U-sväng mot {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Byt till {modifier} körfält",
                "name": "Byt till {modifier} körfält, in på {way_name}",
                "destination": "Byt till {modifier} körfält, mot {destination}"
            },
            "straight": {
                "default": "Fortsätt",
                "name": "Kör in på {way_name}",
                "destination": "Kör mot {destination}"
            },
            "slight left": {
                "default": "Byt till vänstra körfältet",
                "name": "Byt till vänstra körfältet, in på {way_name}",
                "destination": "Byt till vänstra körfältet, mot {destination}"
            },
            "slight right": {
                "default": "Byt till högra körfältet",
                "name": "Byt till högra körfältet, in på {way_name}",
                "destination": "Byt till högra körfältet, mot {destination}"
            },
            "sharp left": {
                "default": "Byt till vänstra körfältet",
                "name": "Byt till vänstra körfältet, in på {way_name}",
                "destination": "Byt till vänstra körfältet, mot {destination}"
            },
            "sharp right": {
                "default": "Byt till högra körfältet",
                "name": "Byt till högra körfältet, in på {way_name}",
                "destination": "Byt till högra körfältet, mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng",
                "name": "Gör en U-sväng in på {way_name}",
                "destination": "Gör en U-sväng mot {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Fortsätt {modifier}",
                "name": "Fortsätt {modifier} på {way_name}",
                "destination": "Fortsätt {modifier} mot {destination}"
            },
            "straight": {
                "default": "Fortsätt rakt fram",
                "name": "Fortsätt in på {way_name}",
                "destination": "Fortsätt mot {destination}"
            },
            "sharp left": {
                "default": "Gör en skarp vänstersväng",
                "name": "Gör en skarp vänstersväng in på {way_name}",
                "destination": "Gör en skarp vänstersväng mot {destination}"
            },
            "sharp right": {
                "default": "Gör en skarp högersväng",
                "name": "Gör en skarp högersväng in på {way_name}",
                "destination": "Gör en skarp högersväng mot {destination}"
            },
            "slight left": {
                "default": "Fortsätt med lätt vänstersväng",
                "name": "Fortsätt med lätt vänstersväng in på {way_name}",
                "destination": "Fortsätt med lätt vänstersväng mot {destination}"
            },
            "slight right": {
                "default": "Fortsätt med lätt högersväng",
                "name": "Fortsätt med lätt högersväng in på {way_name}",
                "destination": "Fortsätt med lätt högersväng mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng",
                "name": "Gör en U-sväng in på {way_name}",
                "destination": "Gör en U-sväng mot {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Fortsätt {modifier}",
                "name": "Fortsätt {modifier} på {way_name}",
                "destination": "Fortsätt {modifier} mot {destination}"
            },
            "uturn": {
                "default": "Gör en U-sväng",
                "name": "Gör en U-sväng in på {way_name}",
                "destination": "Gör en U-sväng mot {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Ta avfarten",
                "name": "Ta avfarten in på {way_name}",
                "destination": "Ta avfarten mot {destination}",
                "exit": "Ta avfart {exit} ",
                "exit_destination": "Ta avfart {exit} mot {destination}"
            },
            "left": {
                "default": "Ta avfarten till vänster",
                "name": "Ta avfarten till vänster in på {way_name}",
                "destination": "Ta avfarten till vänster mot {destination}",
                "exit": "Ta avfart {exit} till vänster",
                "exit_destination": "Ta avfart {exit} till vänster mot {destination}"
            },
            "right": {
                "default": "Ta avfarten till höger",
                "name": "Ta avfarten till höger in på {way_name}",
                "destination": "Ta avfarten till höger mot {destination}",
                "exit": "Ta avfart {exit} till höger",
                "exit_destination": "Ta avfart {exit} till höger mot {destination}"
            },
            "sharp left": {
                "default": "Ta avfarten till vänster",
                "name": "Ta avfarten till vänster in på {way_name}",
                "destination": "Ta avfarten till vänster mot {destination}",
                "exit": "Ta avfart {exit} till vänster",
                "exit_destination": "Ta avfart {exit} till vänster mot {destination}"
            },
            "sharp right": {
                "default": "Ta avfarten till höger",
                "name": "Ta avfarten till höger in på {way_name}",
                "destination": "Ta avfarten till höger mot {destination}",
                "exit": "Ta avfart {exit} till höger",
                "exit_destination": "Ta avfart {exit} till höger mot {destination}"
            },
            "slight left": {
                "default": "Ta avfarten till vänster",
                "name": "Ta avfarten till vänster in på {way_name}",
                "destination": "Ta avfarten till vänster mot {destination}",
                "exit": "Ta avfart {exit} till vänster",
                "exit_destination": "Ta avfart{exit} till vänster mot {destination}"
            },
            "slight right": {
                "default": "Ta avfarten till höger",
                "name": "Ta avfarten till höger in på {way_name}",
                "destination": "Ta avfarten till höger mot {destination}",
                "exit": "Ta avfart {exit} till höger",
                "exit_destination": "Ta avfart {exit} till höger mot {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Ta påfarten",
                "name": "Ta påfarten in på {way_name}",
                "destination": "Ta påfarten mot {destination}"
            },
            "left": {
                "default": "Ta påfarten till vänster",
                "name": "Ta påfarten till vänster in på {way_name}",
                "destination": "Ta påfarten till vänster mot {destination}"
            },
            "right": {
                "default": "Ta påfarten till höger",
                "name": "Ta påfarten till höger in på {way_name}",
                "destination": "Ta påfarten till höger mot {destination}"
            },
            "sharp left": {
                "default": "Ta påfarten till vänster",
                "name": "Ta påfarten till vänster in på {way_name}",
                "destination": "Ta påfarten till vänster mot {destination}"
            },
            "sharp right": {
                "default": "Ta påfarten till höger",
                "name": "Ta påfarten till höger in på {way_name}",
                "destination": "Ta påfarten till höger mot {destination}"
            },
            "slight left": {
                "default": "Ta påfarten till vänster",
                "name": "Ta påfarten till vänster in på {way_name}",
                "destination": "Ta påfarten till vänster mot {destination}"
            },
            "slight right": {
                "default": "Ta påfarten till höger",
                "name": "Ta påfarten till höger in på {way_name}",
                "destination": "Ta påfarten till höger mot {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Kör in i rondellen",
                    "name": "I rondellen, ta avfarten in på {way_name}",
                    "destination": "I rondellen, ta av mot {destination}"
                },
                "name": {
                    "default": "Kör in i {rotary_name}",
                    "name": "I {rotary_name}, ta av in på {way_name}",
                    "destination": "I {rotary_name}, ta av mot {destination}"
                },
                "exit": {
                    "default": "I rondellen, ta {exit_number} avfarten",
                    "name": "I rondellen, ta {exit_number} avfarten in på {way_name}",
                    "destination": "I rondellen, ta {exit_number} avfarten mot {destination}"
                },
                "name_exit": {
                    "default": "I {rotary_name}, ta {exit_number} avfarten",
                    "name": "I {rotary_name}, ta {exit_number}  avfarten in på {way_name}",
                    "destination": "I {rotary_name}, ta {exit_number} avfarten mot {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "I rondellen, ta {exit_number} avfarten",
                    "name": "I rondellen, ta {exit_number} avfarten in på {way_name}",
                    "destination": "I rondellen, ta {exit_number} avfarten mot {destination}"
                },
                "default": {
                    "default": "Kör in i rondellen",
                    "name": "I rondellen, ta avfarten in på {way_name}",
                    "destination": "I rondellen, ta av mot {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Sväng {modifier}",
                "name": "Sväng {modifier} in på {way_name}",
                "destination": "Sväng {modifier} mot {destination}"
            },
            "left": {
                "default": "Sväng vänster",
                "name": "Sväng vänster in på {way_name}",
                "destination": "Sväng vänster mot {destination}"
            },
            "right": {
                "default": "Sväng höger",
                "name": "Sväng höger in på {way_name}",
                "destination": "Sväng höger mot {destination}"
            },
            "straight": {
                "default": "Fortsätt rakt fram",
                "name": "Fortsätt rakt fram in på {way_name}",
                "destination": "Fortsätt rakt fram mot {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Kör ut ur rondellen",
                "name": "Kör ut ur rondellen in på {way_name}",
                "destination": "Kör ut ur rondellen mot {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Kör ut ur rondellen",
                "name": "Kör ut ur rondellen in på {way_name}",
                "destination": "Kör ut ur rondellen mot {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Sväng {modifier}",
                "name": "Sväng {modifier} in på {way_name}",
                "destination": "Sväng {modifier} mot {destination}"
            },
            "left": {
                "default": "Sväng vänster",
                "name": "Sväng vänster in på {way_name}",
                "destination": "Sväng vänster mot {destination}"
            },
            "right": {
                "default": "Sväng höger",
                "name": "Sväng höger in på {way_name}",
                "destination": "Sväng höger mot {destination}"
            },
            "straight": {
                "default": "Kör rakt fram",
                "name": "Kör rakt fram in på {way_name}",
                "destination": "Kör rakt fram mot {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Fortsätt rakt fram"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],44:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "birinci",
                "2": "ikinci",
                "3": "üçüncü",
                "4": "dördüncü",
                "5": "beşinci",
                "6": "altıncı",
                "7": "yedinci",
                "8": "sekizinci",
                "9": "dokuzuncu",
                "10": "onuncu"
            },
            "direction": {
                "north": "kuzey",
                "northeast": "kuzeydoğu",
                "east": "doğu",
                "southeast": "güneydoğu",
                "south": "güney",
                "southwest": "güneybatı",
                "west": "batı",
                "northwest": "kuzeybatı"
            },
            "modifier": {
                "left": "sol",
                "right": "sağ",
                "sharp left": "keskin sol",
                "sharp right": "keskin sağ",
                "slight left": "hafif sol",
                "slight right": "hafif sağ",
                "straight": "düz",
                "uturn": "U dönüşü"
            },
            "lanes": {
                "xo": "Sağda kalın",
                "ox": "Solda kalın",
                "xox": "Ortada kalın",
                "oxo": "Solda veya sağda kalın"
            }
        },
        "modes": {
            "ferry": {
                "default": "Vapur kullan",
                "name": "{way_name} vapurunu kullan",
                "destination": "{destination} istikametine giden vapuru kullan"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one} ve {distance} sonra {instruction_two}",
            "two linked": "{instruction_one} ve sonra {instruction_two}",
            "one in distance": "{distance} sonra, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "exit {exit}"
        },
        "arrive": {
            "default": {
                "default": "{nth} hedefinize ulaştınız",
                "upcoming": "{nth} hedefinize ulaştınız",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız"
            },
            "left": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz solunuzdadır"
            },
            "right": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz sağınızdadır"
            },
            "sharp left": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz solunuzdadır"
            },
            "sharp right": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz sağınızdadır"
            },
            "slight right": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz sağınızdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz sağınızdadır"
            },
            "slight left": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz solunuzdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz solunuzdadır"
            },
            "straight": {
                "default": "{nth} hedefinize ulaştınız, hedefiniz karşınızdadır",
                "upcoming": "{nth} hedefinize ulaştınız, hedefiniz karşınızdadır",
                "short": "{nth} hedefinize ulaştınız",
                "short-upcoming": "{nth} hedefinize ulaştınız",
                "named": "{waypoint_name} ulaştınız, hedefiniz karşınızdadır"
            }
        },
        "continue": {
            "default": {
                "default": "{modifier} yöne dön",
                "name": "{way_name} üzerinde kalmak için {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön",
                "exit": "{way_name} üzerinde {modifier} yöne dön"
            },
            "straight": {
                "default": "Düz devam edin",
                "name": "{way_name} üzerinde kalmak için düz devam et",
                "destination": "{destination} istikametinde devam et",
                "distance": "{distance} boyunca düz devam et",
                "namedistance": "{distance} boyunca {way_name} üzerinde devam et"
            },
            "sharp left": {
                "default": "Sola keskin dönüş yap",
                "name": "{way_name} üzerinde kalmak için sola keskin dönüş yap",
                "destination": "{destination} istikametinde sola keskin dönüş yap"
            },
            "sharp right": {
                "default": "Sağa keskin dönüş yap",
                "name": "{way_name} üzerinde kalmak için sağa keskin dönüş yap",
                "destination": "{destination} istikametinde sağa keskin dönüş yap"
            },
            "slight left": {
                "default": "Sola hafif dönüş yap",
                "name": "{way_name} üzerinde kalmak için sola hafif dönüş yap",
                "destination": "{destination} istikametinde sola hafif dönüş yap"
            },
            "slight right": {
                "default": "Sağa hafif dönüş yap",
                "name": "{way_name} üzerinde kalmak için sağa hafif dönüş yap",
                "destination": "{destination} istikametinde sağa hafif dönüş yap"
            },
            "uturn": {
                "default": "U dönüşü yapın",
                "name": "Bir U-dönüşü yap ve {way_name} devam et",
                "destination": "{destination} istikametinde bir U-dönüşü yap"
            }
        },
        "depart": {
            "default": {
                "default": "{direction} tarafına yönelin",
                "name": "{way_name} üzerinde {direction} yöne git",
                "namedistance": "Head {direction} on {way_name} for {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "{modifier} tarafa dönün",
                "name": "{way_name} üzerinde {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön"
            },
            "straight": {
                "default": "Düz devam edin",
                "name": "{way_name} üzerinde düz devam et",
                "destination": "{destination} istikametinde düz devam et"
            },
            "uturn": {
                "default": "Yolun sonunda U dönüşü yapın",
                "name": "Yolun sonunda {way_name} üzerinde bir U-dönüşü yap",
                "destination": "Yolun sonunda {destination} istikametinde bir U-dönüşü yap"
            }
        },
        "fork": {
            "default": {
                "default": "Yol ayrımında {modifier} yönde kal",
                "name": "{way_name} üzerindeki yol ayrımında {modifier} yönde kal",
                "destination": "{destination} istikametindeki yol ayrımında {modifier} yönde kal"
            },
            "slight left": {
                "default": "Çatalın solundan devam edin",
                "name": "Çatalın solundan {way_name} yoluna doğru ",
                "destination": "{destination} istikametindeki yol ayrımında solda kal"
            },
            "slight right": {
                "default": "Çatalın sağından devam edin",
                "name": "{way_name} üzerindeki yol ayrımında sağda kal",
                "destination": "{destination} istikametindeki yol ayrımında sağda kal"
            },
            "sharp left": {
                "default": "Çatalda keskin sola dönün",
                "name": "{way_name} yoluna doğru sola keskin dönüş yapın",
                "destination": "{destination} istikametinde sola keskin dönüş yap"
            },
            "sharp right": {
                "default": "Çatalda keskin sağa dönün",
                "name": "{way_name} yoluna doğru sağa keskin dönüş yapın",
                "destination": "{destination} istikametinde sağa keskin dönüş yap"
            },
            "uturn": {
                "default": "U dönüşü yapın",
                "name": "{way_name} yoluna U dönüşü yapın",
                "destination": "{destination} istikametinde bir U-dönüşü yap"
            }
        },
        "merge": {
            "default": {
                "default": "{modifier} yöne gir",
                "name": "{way_name} üzerinde {modifier} yöne gir",
                "destination": "{destination} istikametinde {modifier} yöne gir"
            },
            "straight": {
                "default": "düz yöne gir",
                "name": "{way_name} üzerinde düz yöne gir",
                "destination": "{destination} istikametinde düz yöne gir"
            },
            "slight left": {
                "default": "Sola gir",
                "name": "{way_name} üzerinde sola gir",
                "destination": "{destination} istikametinde sola gir"
            },
            "slight right": {
                "default": "Sağa gir",
                "name": "{way_name} üzerinde sağa gir",
                "destination": "{destination} istikametinde sağa gir"
            },
            "sharp left": {
                "default": "Sola gir",
                "name": "{way_name} üzerinde sola gir",
                "destination": "{destination} istikametinde sola gir"
            },
            "sharp right": {
                "default": "Sağa gir",
                "name": "{way_name} üzerinde sağa gir",
                "destination": "{destination} istikametinde sağa gir"
            },
            "uturn": {
                "default": "U dönüşü yapın",
                "name": "{way_name} yoluna U dönüşü yapın",
                "destination": "{destination} istikametinde bir U-dönüşü yap"
            }
        },
        "new name": {
            "default": {
                "default": "{modifier} yönde devam et",
                "name": "{way_name} üzerinde {modifier} yönde devam et",
                "destination": "{destination} istikametinde {modifier} yönde devam et"
            },
            "straight": {
                "default": "Düz devam et",
                "name": "{way_name} üzerinde devam et",
                "destination": "{destination} istikametinde devam et"
            },
            "sharp left": {
                "default": "Sola keskin dönüş yapın",
                "name": "{way_name} yoluna doğru sola keskin dönüş yapın",
                "destination": "{destination} istikametinde sola keskin dönüş yap"
            },
            "sharp right": {
                "default": "Sağa keskin dönüş yapın",
                "name": "{way_name} yoluna doğru sağa keskin dönüş yapın",
                "destination": "{destination} istikametinde sağa keskin dönüş yap"
            },
            "slight left": {
                "default": "Hafif soldan devam edin",
                "name": "{way_name} üzerinde hafif solda devam et",
                "destination": "{destination} istikametinde hafif solda devam et"
            },
            "slight right": {
                "default": "Hafif sağdan devam edin",
                "name": "{way_name} üzerinde hafif sağda devam et",
                "destination": "{destination} istikametinde hafif sağda devam et"
            },
            "uturn": {
                "default": "U dönüşü yapın",
                "name": "{way_name} yoluna U dönüşü yapın",
                "destination": "{destination} istikametinde bir U-dönüşü yap"
            }
        },
        "notification": {
            "default": {
                "default": "{modifier} yönde devam et",
                "name": "{way_name} üzerinde {modifier} yönde devam et",
                "destination": "{destination} istikametinde {modifier} yönde devam et"
            },
            "uturn": {
                "default": "U dönüşü yapın",
                "name": "{way_name} yoluna U dönüşü yapın",
                "destination": "{destination} istikametinde bir U-dönüşü yap"
            }
        },
        "off ramp": {
            "default": {
                "default": "Bağlantı yoluna geç",
                "name": "{way_name} üzerindeki bağlantı yoluna geç",
                "destination": "{destination} istikametine giden bağlantı yoluna geç",
                "exit": "{exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} çıkış yoluna geç"
            },
            "left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç",
                "exit": "Soldaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sol çıkış yoluna geç"
            },
            "right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç",
                "exit": "Sağdaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sağ çıkış yoluna geç"
            },
            "sharp left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç",
                "exit": "Soldaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sol çıkış yoluna geç"
            },
            "sharp right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç",
                "exit": "Sağdaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sağ çıkış yoluna geç"
            },
            "slight left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç",
                "exit": "Soldaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sol çıkış yoluna geç"
            },
            "slight right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç",
                "exit": "Sağdaki {exit} çıkış yoluna geç",
                "exit_destination": "{destination} istikametindeki {exit} sağ çıkış yoluna geç"
            }
        },
        "on ramp": {
            "default": {
                "default": "Bağlantı yoluna geç",
                "name": "{way_name} üzerindeki bağlantı yoluna geç",
                "destination": "{destination} istikametine giden bağlantı yoluna geç"
            },
            "left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç"
            },
            "right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç"
            },
            "sharp left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç"
            },
            "sharp right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç"
            },
            "slight left": {
                "default": "Soldaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sol bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sol bağlantı yoluna geç"
            },
            "slight right": {
                "default": "Sağdaki bağlantı yoluna geç",
                "name": "{way_name} üzerindeki sağ bağlantı yoluna geç",
                "destination": "{destination} istikametine giden sağ bağlantı yoluna geç"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Dönel kavşağa gir",
                    "name": "Dönel kavşağa gir ve {way_name} üzerinde çık",
                    "destination": "Dönel kavşağa gir ve {destination} istikametinde çık"
                },
                "name": {
                    "default": "{rotary_name} dönel kavşağa gir",
                    "name": "{rotary_name} dönel kavşağa gir ve {way_name} üzerinde çık",
                    "destination": "{rotary_name} dönel kavşağa gir ve {destination} istikametinde çık"
                },
                "exit": {
                    "default": "Dönel kavşağa gir ve {exit_number} numaralı çıkışa gir",
                    "name": "Dönel kavşağa gir ve {way_name} üzerindeki {exit_number} numaralı çıkışa gir",
                    "destination": "Dönel kavşağa gir ve {destination} istikametindeki {exit_number} numaralı çıkışa gir"
                },
                "name_exit": {
                    "default": "{rotary_name} dönel kavşağa gir ve {exit_number} numaralı çıkışa gir",
                    "name": "{rotary_name} dönel kavşağa gir ve {way_name} üzerindeki {exit_number} numaralı çıkışa gir",
                    "destination": "{rotary_name} dönel kavşağa gir ve {destination} istikametindeki {exit_number} numaralı çıkışa gir"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Göbekli kavşağa gir ve {exit_number} numaralı çıkışa gir",
                    "name": "Göbekli kavşağa gir ve {way_name} üzerindeki {exit_number} numaralı çıkışa gir",
                    "destination": "Göbekli kavşağa gir ve {destination} istikametindeki {exit_number} numaralı çıkışa gir"
                },
                "default": {
                    "default": "Göbekli kavşağa gir",
                    "name": "Göbekli kavşağa gir ve {way_name} üzerinde çık",
                    "destination": "Göbekli kavşağa gir ve {destination} istikametinde çık"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier} yöne dön",
                "name": "{way_name} üzerinde {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön"
            },
            "left": {
                "default": "Sola dön",
                "name": "{way_name} üzerinde sola dön",
                "destination": "{destination} istikametinde sola dön"
            },
            "right": {
                "default": "Sağa dön",
                "name": "{way_name} üzerinde sağa dön",
                "destination": "{destination} istikametinde sağa dön"
            },
            "straight": {
                "default": "Düz devam et",
                "name": "{way_name} üzerinde düz devam et",
                "destination": "{destination} istikametinde düz devam et"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "{modifier} yöne dön",
                "name": "{way_name} üzerinde {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön"
            },
            "left": {
                "default": "Sola dön",
                "name": "{way_name} üzerinde sola dön",
                "destination": "{destination} istikametinde sola dön"
            },
            "right": {
                "default": "Sağa dön",
                "name": "{way_name} üzerinde sağa dön",
                "destination": "{destination} istikametinde sağa dön"
            },
            "straight": {
                "default": "Düz devam et",
                "name": "{way_name} üzerinde düz devam et",
                "destination": "{destination} istikametinde düz devam et"
            }
        },
        "exit rotary": {
            "default": {
                "default": "{modifier} yöne dön",
                "name": "{way_name} üzerinde {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön"
            },
            "left": {
                "default": "Sola dön",
                "name": "{way_name} üzerinde sola dön",
                "destination": "{destination} istikametinde sola dön"
            },
            "right": {
                "default": "Sağa dön",
                "name": "{way_name} üzerinde sağa dön",
                "destination": "{destination} istikametinde sağa dön"
            },
            "straight": {
                "default": "Düz devam et",
                "name": "{way_name} üzerinde düz devam et",
                "destination": "{destination} istikametinde düz devam et"
            }
        },
        "turn": {
            "default": {
                "default": "{modifier} yöne dön",
                "name": "{way_name} üzerinde {modifier} yöne dön",
                "destination": "{destination} istikametinde {modifier} yöne dön"
            },
            "left": {
                "default": "Sola dönün",
                "name": "{way_name} üzerinde sola dön",
                "destination": "{destination} istikametinde sola dön"
            },
            "right": {
                "default": "Sağa dönün",
                "name": "{way_name} üzerinde sağa dön",
                "destination": "{destination} istikametinde sağa dön"
            },
            "straight": {
                "default": "Düz git",
                "name": "{way_name} üzerinde düz git",
                "destination": "{destination} istikametinde düz git"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Düz devam edin"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],45:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "1й",
                "2": "2й",
                "3": "3й",
                "4": "4й",
                "5": "5й",
                "6": "6й",
                "7": "7й",
                "8": "8й",
                "9": "9й",
                "10": "10й"
            },
            "direction": {
                "north": "північ",
                "northeast": "північний схід",
                "east": "схід",
                "southeast": "південний схід",
                "south": "південь",
                "southwest": "південний захід",
                "west": "захід",
                "northwest": "північний захід"
            },
            "modifier": {
                "left": "ліворуч",
                "right": "праворуч",
                "sharp left": "різко ліворуч",
                "sharp right": "різко праворуч",
                "slight left": "плавно ліворуч",
                "slight right": "плавно праворуч",
                "straight": "прямо",
                "uturn": "розворот"
            },
            "lanes": {
                "xo": "Тримайтесь праворуч",
                "ox": "Тримайтесь ліворуч",
                "xox": "Тримайтесь в середині",
                "oxo": "Тримайтесь праворуч або ліворуч"
            }
        },
        "modes": {
            "ferry": {
                "default": "Скористайтесь поромом",
                "name": "Скористайтесь поромом {way_name}",
                "destination": "Скористайтесь поромом у напрямку {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, потім, через {distance}, {instruction_two}",
            "two linked": "{instruction_one}, потім {instruction_two}",
            "one in distance": "Через {distance}, {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "з'їзд {exit}"
        },
        "arrive": {
            "default": {
                "default": "Ви прибули у ваш {nth} пункт призначення",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name}"
            },
            "left": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – ліворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, ліворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} ліворуч"
            },
            "right": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – праворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, праворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} праворуч"
            },
            "sharp left": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – ліворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, ліворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} ліворуч"
            },
            "sharp right": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – праворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, праворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} праворуч"
            },
            "slight right": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – праворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, праворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} праворуч"
            },
            "slight left": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – ліворуч",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, ліворуч",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} ліворуч"
            },
            "straight": {
                "default": "Ви прибули у ваш {nth} пункт призначення, він – прямо перед вами",
                "upcoming": "Ви наближаєтесь до вашого {nth} місця призначення, прямо перед вами",
                "short": "Ви прибули",
                "short-upcoming": "Ви прибудете",
                "named": "Ви прибули у {waypoint_name} прямо перед вами"
            }
        },
        "continue": {
            "default": {
                "default": "Поверніть {modifier}",
                "name": "Поверніть{modifier} залишаючись на {way_name}",
                "destination": "Поверніть {modifier} у напрямку {destination}",
                "exit": "Поверніть {modifier} на {way_name}"
            },
            "straight": {
                "default": "Продовжуйте рух прямо",
                "name": "Продовжуйте рух прямо залишаючись на {way_name}",
                "destination": "Рухайтесь у напрямку {destination}",
                "distance": "Продовжуйте рух прямо {distance}",
                "namedistance": "Продовжуйте рух по {way_name} {distance}"
            },
            "sharp left": {
                "default": "Поверніть різко ліворуч",
                "name": "Поверніть різко ліворуч щоб залишитись на {way_name}",
                "destination": "Поверніть різко ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "Поверніть різко праворуч",
                "name": "Поверніть різко праворуч щоб залишитись на {way_name}",
                "destination": "Поверніть різко праворуч у напрямку {destination}"
            },
            "slight left": {
                "default": "Поверніть різко ліворуч",
                "name": "Поверніть плавно ліворуч щоб залишитись на {way_name}",
                "destination": "Поверніть плавно ліворуч у напрямку {destination}"
            },
            "slight right": {
                "default": "Поверніть плавно праворуч",
                "name": "Поверніть плавно праворуч щоб залишитись на {way_name}",
                "destination": "Поверніть плавно праворуч у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот",
                "name": "Здійсніть розворот та рухайтесь по {way_name}",
                "destination": "Здійсніть розворот у напрямку {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Прямуйте на {direction}",
                "name": "Прямуйте на {direction} по {way_name}",
                "namedistance": "Прямуйте на {direction} по {way_name} {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Поверніть {modifier}",
                "name": "Поверніть {modifier} на {way_name}",
                "destination": "Поверніть {modifier} у напрямку {destination}"
            },
            "straight": {
                "default": "Продовжуйте рух прямо",
                "name": "Продовжуйте рух прямо до {way_name}",
                "destination": "Продовжуйте рух прямо у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот в кінці дороги",
                "name": "Здійсніть розворот на {way_name} в кінці дороги",
                "destination": "Здійсніть розворот у напрямку {destination} в кінці дороги"
            }
        },
        "fork": {
            "default": {
                "default": "На роздоріжжі тримайтеся {modifier}",
                "name": "Тримайтеся {modifier} і рухайтесь на {way_name}",
                "destination": "Тримайтеся {modifier} в напрямку {destination}"
            },
            "slight left": {
                "default": "На роздоріжжі тримайтеся ліворуч",
                "name": "Тримайтеся ліворуч і рухайтесь на {way_name}",
                "destination": "Тримайтеся ліворуч в напрямку {destination}"
            },
            "slight right": {
                "default": "На роздоріжжі тримайтеся праворуч",
                "name": "Тримайтеся праворуч і рухайтесь на {way_name}",
                "destination": "Тримайтеся праворуч в напрямку {destination}"
            },
            "sharp left": {
                "default": "На роздоріжжі різко поверніть ліворуч",
                "name": "Прийміть різко ліворуч на {way_name}",
                "destination": "Прийміть різко ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "На роздоріжжі різко поверніть праворуч",
                "name": "Прийміть різко праворуч на {way_name}",
                "destination": "Прийміть різко праворуч у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот",
                "name": "Здійсніть розворот на {way_name}",
                "destination": "Здійсніть розворот у напрямку {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Приєднайтеся до потоку {modifier}",
                "name": "Приєднайтеся до потоку {modifier} на {way_name}",
                "destination": "Приєднайтеся до потоку {modifier} у напрямку {destination}"
            },
            "straight": {
                "default": "Приєднайтеся до потоку",
                "name": "Приєднайтеся до потоку на {way_name}",
                "destination": "Приєднайтеся до потоку у напрямку {destination}"
            },
            "slight left": {
                "default": "Приєднайтеся до потоку ліворуч",
                "name": "Приєднайтеся до потоку ліворуч на {way_name}",
                "destination": "Приєднайтеся до потоку ліворуч у напрямку {destination}"
            },
            "slight right": {
                "default": "Приєднайтеся до потоку праворуч",
                "name": "Приєднайтеся до потоку праворуч на {way_name}",
                "destination": "Приєднайтеся до потоку праворуч у напрямку {destination}"
            },
            "sharp left": {
                "default": "Приєднайтеся до потоку ліворуч",
                "name": "Приєднайтеся до потоку ліворуч на {way_name}",
                "destination": "Приєднайтеся до потоку ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "Приєднайтеся до потоку праворуч",
                "name": "Приєднайтеся до потоку праворуч на {way_name}",
                "destination": "Приєднайтеся до потоку праворуч у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот",
                "name": "Здійсніть розворот на {way_name}",
                "destination": "Здійсніть розворот у напрямку {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Рухайтесь {modifier}",
                "name": "Рухайтесь {modifier} на {way_name}",
                "destination": "Рухайтесь {modifier} у напрямку {destination}"
            },
            "straight": {
                "default": "Рухайтесь прямо",
                "name": "Рухайтесь по {way_name}",
                "destination": "Рухайтесь у напрямку {destination}"
            },
            "sharp left": {
                "default": "Прийміть різко ліворуч",
                "name": "Прийміть різко ліворуч на {way_name}",
                "destination": "Прийміть різко ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "Прийміть різко праворуч",
                "name": "Прийміть різко праворуч на {way_name}",
                "destination": "Прийміть різко праворуч у напрямку {destination}"
            },
            "slight left": {
                "default": "Рухайтесь плавно ліворуч",
                "name": "Рухайтесь плавно ліворуч на {way_name}",
                "destination": "Рухайтесь плавно ліворуч у напрямку {destination}"
            },
            "slight right": {
                "default": "Рухайтесь плавно праворуч",
                "name": "Рухайтесь плавно праворуч на {way_name}",
                "destination": "Рухайтесь плавно праворуч у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот",
                "name": "Здійсніть розворот на {way_name}",
                "destination": "Здійсніть розворот у напрямку {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Рухайтесь {modifier}",
                "name": "Рухайтесь {modifier} на {way_name}",
                "destination": "Рухайтесь {modifier} у напрямку {destination}"
            },
            "uturn": {
                "default": "Здійсніть розворот",
                "name": "Здійсніть розворот на {way_name}",
                "destination": "Здійсніть розворот у напрямку {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Рухайтесь на зʼїзд",
                "name": "Рухайтесь на зʼїзд на {way_name}",
                "destination": "Рухайтесь на зʼїзд у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit}",
                "exit_destination": "Оберіть з'їзд {exit} у напрямку {destination}"
            },
            "left": {
                "default": "Рухайтесь на зʼїзд ліворуч",
                "name": "Рухайтесь на зʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд ліворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} ліворуч",
                "exit_destination": "Оберіть з'їзд {exit} ліворуч у напрямку {destination}"
            },
            "right": {
                "default": "Рухайтесь на зʼїзд праворуч",
                "name": "Рухайтесь на зʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд праворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} праворуч",
                "exit_destination": "Оберіть з'їзд {exit} праворуч у напрямку {destination}"
            },
            "sharp left": {
                "default": "Рухайтесь на зʼїзд ліворуч",
                "name": "Рухайтесь на зʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд ліворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} ліворуч",
                "exit_destination": "Оберіть з'їзд {exit} ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "Рухайтесь на зʼїзд праворуч",
                "name": "Рухайтесь на зʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд праворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} праворуч",
                "exit_destination": "Оберіть з'їзд {exit} праворуч у напрямку {destination}"
            },
            "slight left": {
                "default": "Рухайтесь на зʼїзд ліворуч",
                "name": "Рухайтесь на зʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд ліворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} ліворуч",
                "exit_destination": "Оберіть з'їзд {exit} ліворуч у напрямку {destination}"
            },
            "slight right": {
                "default": "Рухайтесь на зʼїзд праворуч",
                "name": "Рухайтесь на зʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на зʼїзд праворуч у напрямку {destination}",
                "exit": "Оберіть з'їзд {exit} праворуч",
                "exit_destination": "Оберіть з'їзд {exit} праворуч у напрямку {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Рухайтесь на вʼїзд",
                "name": "Рухайтесь на вʼїзд на {way_name}",
                "destination": "Рухайтесь на вʼїзд у напрямку {destination}"
            },
            "left": {
                "default": "Рухайтесь на вʼїзд ліворуч",
                "name": "Рухайтесь на вʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд ліворуч у напрямку {destination}"
            },
            "right": {
                "default": "Рухайтесь на вʼїзд праворуч",
                "name": "Рухайтесь на вʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд праворуч у напрямку {destination}"
            },
            "sharp left": {
                "default": "Рухайтесь на вʼїзд ліворуч",
                "name": "Рухайтесь на вʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд ліворуч у напрямку {destination}"
            },
            "sharp right": {
                "default": "Рухайтесь на вʼїзд праворуч",
                "name": "Рухайтесь на вʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд праворуч у напрямку {destination}"
            },
            "slight left": {
                "default": "Рухайтесь на вʼїзд ліворуч",
                "name": "Рухайтесь на вʼїзд ліворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд ліворуч у напрямку {destination}"
            },
            "slight right": {
                "default": "Рухайтесь на вʼїзд праворуч",
                "name": "Рухайтесь на вʼїзд праворуч на {way_name}",
                "destination": "Рухайтесь на вʼїзд праворуч у напрямку {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Рухайтесь по колу",
                    "name": "Рухайтесь по колу до {way_name}",
                    "destination": "Рухайтесь по колу в напрямку {destination}"
                },
                "name": {
                    "default": "Рухайтесь по {rotary_name}",
                    "name": "Рухайтесь по {rotary_name} та поверніть на {way_name}",
                    "destination": "Рухайтесь по {rotary_name} та поверніть в напрямку {destination}"
                },
                "exit": {
                    "default": "Рухайтесь по колу та повереніть у {exit_number} з'їзд",
                    "name": "Рухайтесь по колу та поверніть у {exit_number} з'їзд на {way_name}",
                    "destination": "Рухайтесь по колу та поверніть у {exit_number} з'їзд у напрямку {destination}"
                },
                "name_exit": {
                    "default": "Рухайтесь по {rotary_name} та поверніть у {exit_number} з'їзд",
                    "name": "Рухайтесь по {rotary_name} та поверніть у {exit_number} з'їзд на {way_name}",
                    "destination": "Рухайтесь по {rotary_name} та поверніть у {exit_number} з'їзд в напрямку {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Рухайтесь по колу та повереніть у {exit_number} з'їзд",
                    "name": "Рухайтесь по колу та поверніть у {exit_number} з'їзд на {way_name}",
                    "destination": "Рухайтесь по колу та поверніть у {exit_number} з'їзд у напрямку {destination}"
                },
                "default": {
                    "default": "Рухайтесь по колу",
                    "name": "Рухайтесь по колу до {way_name}",
                    "destination": "Рухайтесь по колу в напрямку {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Рухайтесь {modifier}",
                "name": "Рухайтесь {modifier} на {way_name}",
                "destination": "Рухайтесь {modifier} в напрямку {destination}"
            },
            "left": {
                "default": "Поверніть ліворуч",
                "name": "Поверніть ліворуч на {way_name}",
                "destination": "Поверніть ліворуч у напрямку {destination}"
            },
            "right": {
                "default": "Поверніть праворуч",
                "name": "Поверніть праворуч на {way_name}",
                "destination": "Поверніть праворуч у напрямку {destination}"
            },
            "straight": {
                "default": "Рухайтесь прямо",
                "name": "Продовжуйте рух прямо до {way_name}",
                "destination": "Продовжуйте рух прямо у напрямку {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Залишить коло",
                "name": "Залишить коло на {way_name} зʼїзді",
                "destination": "Залишить коло в напрямку {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Залишить коло",
                "name": "Залишить коло на {way_name} зʼїзді",
                "destination": "Залишить коло в напрямку {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Рухайтесь {modifier}",
                "name": "Рухайтесь {modifier} на {way_name}",
                "destination": "Рухайтесь {modifier} в напрямку {destination}"
            },
            "left": {
                "default": "Поверніть ліворуч",
                "name": "Поверніть ліворуч на {way_name}",
                "destination": "Поверніть ліворуч у напрямку {destination}"
            },
            "right": {
                "default": "Поверніть праворуч",
                "name": "Поверніть праворуч на {way_name}",
                "destination": "Поверніть праворуч у напрямку {destination}"
            },
            "straight": {
                "default": "Рухайтесь прямо",
                "name": "Рухайтесь прямо по {way_name}",
                "destination": "Рухайтесь прямо у напрямку {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Продовжуйте рух прямо"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],46:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": true
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "đầu tiên",
                "2": "thứ 2",
                "3": "thứ 3",
                "4": "thứ 4",
                "5": "thứ 5",
                "6": "thú 6",
                "7": "thứ 7",
                "8": "thứ 8",
                "9": "thứ 9",
                "10": "thứ 10"
            },
            "direction": {
                "north": "bắc",
                "northeast": "đông bắc",
                "east": "đông",
                "southeast": "đông nam",
                "south": "nam",
                "southwest": "tây nam",
                "west": "tây",
                "northwest": "tây bắc"
            },
            "modifier": {
                "left": "trái",
                "right": "phải",
                "sharp left": "trái gắt",
                "sharp right": "phải gắt",
                "slight left": "trái nghiêng",
                "slight right": "phải nghiêng",
                "straight": "thẳng",
                "uturn": "ngược"
            },
            "lanes": {
                "xo": "Đi bên phải",
                "ox": "Đi bên trái",
                "xox": "Đi vào giữa",
                "oxo": "Đi bên trái hay bên phải"
            }
        },
        "modes": {
            "ferry": {
                "default": "Lên phà",
                "name": "Lên phà {way_name}",
                "destination": "Lên phà đi {destination}"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}, rồi {distance} nữa thì {instruction_two}",
            "two linked": "{instruction_one}, rồi {instruction_two}",
            "one in distance": "{distance} nữa thì {instruction_one}",
            "name and ref": "{name} ({ref})",
            "exit with number": "lối ra {exit}"
        },
        "arrive": {
            "default": {
                "default": "Đến nơi {nth}",
                "upcoming": "Đến nơi {nth}",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name}"
            },
            "left": {
                "default": "Đến nơi {nth} ở bên trái",
                "upcoming": "Đến nơi {nth} ở bên trái",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên trái"
            },
            "right": {
                "default": "Đến nơi {nth} ở bên phải",
                "upcoming": "Đến nơi {nth} ở bên phải",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên phải"
            },
            "sharp left": {
                "default": "Đến nơi {nth} ở bên trái",
                "upcoming": "Đến nơi {nth} ở bên trái",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên trái"
            },
            "sharp right": {
                "default": "Đến nơi {nth} ở bên phải",
                "upcoming": "Đến nơi {nth} ở bên phải",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên phải"
            },
            "slight right": {
                "default": "Đến nơi {nth} ở bên phải",
                "upcoming": "Đến nơi {nth} ở bên phải",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên phải"
            },
            "slight left": {
                "default": "Đến nơi {nth} ở bên trái",
                "upcoming": "Đến nơi {nth} ở bên trái",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở bên trái"
            },
            "straight": {
                "default": "Đến nơi {nth} ở trước mặt",
                "upcoming": "Đến nơi {nth} ở trước mặt",
                "short": "Đến nơi",
                "short-upcoming": "Đến nơi",
                "named": "Đến {waypoint_name} ở trước mặt"
            }
        },
        "continue": {
            "default": {
                "default": "Quẹo {modifier}",
                "name": "Quẹo {modifier} để chạy tiếp trên {way_name}",
                "destination": "Quẹo {modifier} về {destination}",
                "exit": "Quẹo {modifier} vào {way_name}"
            },
            "straight": {
                "default": "Chạy thẳng",
                "name": "Chạy tiếp trên {way_name}",
                "destination": "Chạy tiếp về {destination}",
                "distance": "Chạy thẳng cho {distance}",
                "namedistance": "Chạy tiếp trên {way_name} cho {distance}"
            },
            "sharp left": {
                "default": "Quẹo gắt bên trái",
                "name": "Quẹo gắt bên trái để chạy tiếp trên {way_name}",
                "destination": "Quẹo gắt bên trái về {destination}"
            },
            "sharp right": {
                "default": "Quẹo gắt bên phải",
                "name": "Quẹo gắt bên phải để chạy tiếp trên {way_name}",
                "destination": "Quẹo gắt bên phải về {destination}"
            },
            "slight left": {
                "default": "Nghiêng về bên trái",
                "name": "Nghiêng về bên trái để chạy tiếp trên {way_name}",
                "destination": "Nghiêng về bên trái về {destination}"
            },
            "slight right": {
                "default": "Nghiêng về bên phải",
                "name": "Nghiêng về bên phải để chạy tiếp trên {way_name}",
                "destination": "Nghiêng về bên phải về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại",
                "name": "Quẹo ngược lại trên {way_name}",
                "destination": "Quẹo ngược về {destination}"
            }
        },
        "depart": {
            "default": {
                "default": "Đi về hướng {direction}",
                "name": "Đi về hướng {direction} trên {way_name}",
                "namedistance": "Đi về hướng {direction} trên {way_name} cho {distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "Quẹo {modifier}",
                "name": "Quẹo {modifier} vào {way_name}",
                "destination": "Quẹo {modifier} về {destination}"
            },
            "straight": {
                "default": "Chạy thẳng",
                "name": "Chạy tiếp trên {way_name}",
                "destination": "Chạy tiếp về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại tại cuối đường",
                "name": "Quẹo ngược vào {way_name} tại cuối đường",
                "destination": "Quẹo ngược về {destination} tại cuối đường"
            }
        },
        "fork": {
            "default": {
                "default": "Đi bên {modifier} ở ngã ba",
                "name": "Giữ bên {modifier} vào {way_name}",
                "destination": "Giữ bên {modifier} về {destination}"
            },
            "slight left": {
                "default": "Nghiêng về bên trái ở ngã ba",
                "name": "Giữ bên trái vào {way_name}",
                "destination": "Giữ bên trái về {destination}"
            },
            "slight right": {
                "default": "Nghiêng về bên phải ở ngã ba",
                "name": "Giữ bên phải vào {way_name}",
                "destination": "Giữ bên phải về {destination}"
            },
            "sharp left": {
                "default": "Quẹo gắt bên trái ở ngã ba",
                "name": "Quẹo gắt bên trái vào {way_name}",
                "destination": "Quẹo gắt bên trái về {destination}"
            },
            "sharp right": {
                "default": "Quẹo gắt bên phải ở ngã ba",
                "name": "Quẹo gắt bên phải vào {way_name}",
                "destination": "Quẹo gắt bên phải về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại",
                "name": "Quẹo ngược lại {way_name}",
                "destination": "Quẹo ngược lại về {destination}"
            }
        },
        "merge": {
            "default": {
                "default": "Nhập sang {modifier}",
                "name": "Nhập sang {modifier} vào {way_name}",
                "destination": "Nhập sang {modifier} về {destination}"
            },
            "straight": {
                "default": "Nhập đường",
                "name": "Nhập vào {way_name}",
                "destination": "Nhập đường về {destination}"
            },
            "slight left": {
                "default": "Nhập sang trái",
                "name": "Nhập sang trái vào {way_name}",
                "destination": "Nhập sang trái về {destination}"
            },
            "slight right": {
                "default": "Nhập sang phải",
                "name": "Nhập sang phải vào {way_name}",
                "destination": "Nhập sang phải về {destination}"
            },
            "sharp left": {
                "default": "Nhập sang trái",
                "name": "Nhập sang trái vào {way_name}",
                "destination": "Nhập sang trái về {destination}"
            },
            "sharp right": {
                "default": "Nhập sang phải",
                "name": "Nhập sang phải vào {way_name}",
                "destination": "Nhập sang phải về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại",
                "name": "Quẹo ngược lại {way_name}",
                "destination": "Quẹo ngược lại về {destination}"
            }
        },
        "new name": {
            "default": {
                "default": "Chạy tiếp bên {modifier}",
                "name": "Chạy tiếp bên {modifier} trên {way_name}",
                "destination": "Chạy tiếp bên {modifier} về {destination}"
            },
            "straight": {
                "default": "Chạy thẳng",
                "name": "Chạy tiếp trên {way_name}",
                "destination": "Chạy tiếp về {destination}"
            },
            "sharp left": {
                "default": "Quẹo gắt bên trái",
                "name": "Quẹo gắt bên trái vào {way_name}",
                "destination": "Quẹo gắt bên trái về {destination}"
            },
            "sharp right": {
                "default": "Quẹo gắt bên phải",
                "name": "Quẹo gắt bên phải vào {way_name}",
                "destination": "Quẹo gắt bên phải về {destination}"
            },
            "slight left": {
                "default": "Nghiêng về bên trái",
                "name": "Nghiêng về bên trái vào {way_name}",
                "destination": "Nghiêng về bên trái về {destination}"
            },
            "slight right": {
                "default": "Nghiêng về bên phải",
                "name": "Nghiêng về bên phải vào {way_name}",
                "destination": "Nghiêng về bên phải về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại",
                "name": "Quẹo ngược lại {way_name}",
                "destination": "Quẹo ngược lại về {destination}"
            }
        },
        "notification": {
            "default": {
                "default": "Chạy tiếp bên {modifier}",
                "name": "Chạy tiếp bên {modifier} trên {way_name}",
                "destination": "Chạy tiếp bên {modifier} về {destination}"
            },
            "uturn": {
                "default": "Quẹo ngược lại",
                "name": "Quẹo ngược lại {way_name}",
                "destination": "Quẹo ngược lại về {destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "Đi đường nhánh",
                "name": "Đi đường nhánh {way_name}",
                "destination": "Đi đường nhánh về {destination}",
                "exit": "Đi theo lối ra {exit}",
                "exit_destination": "Đi theo lối ra {exit} về {destination}"
            },
            "left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}",
                "exit": "Đi theo lối ra {exit} bên trái",
                "exit_destination": "Đi theo lối ra {exit} bên trái về {destination}"
            },
            "right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}",
                "exit": "Đi theo lối ra {exit} bên phải",
                "exit_destination": "Đi theo lối ra {exit} bên phải về {destination}"
            },
            "sharp left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}",
                "exit": "Đi theo lối ra {exit} bên trái",
                "exit_destination": "Đi theo lối ra {exit} bên trái về {destination}"
            },
            "sharp right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}",
                "exit": "Đi theo lối ra {exit} bên phải",
                "exit_destination": "Đi theo lối ra {exit} bên phải về {destination}"
            },
            "slight left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}",
                "exit": "Đi theo lối ra {exit} bên trái",
                "exit_destination": "Đi theo lối ra {exit} bên trái về {destination}"
            },
            "slight right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}",
                "exit": "Đi theo lối ra {exit} bên phải",
                "exit_destination": "Đi theo lối ra {exit} bên phải về {destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "Đi đường nhánh",
                "name": "Đi đường nhánh {way_name}",
                "destination": "Đi đường nhánh về {destination}"
            },
            "left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}"
            },
            "right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}"
            },
            "sharp left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}"
            },
            "sharp right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}"
            },
            "slight left": {
                "default": "Đi đường nhánh bên trái",
                "name": "Đi đường nhánh {way_name} bên trái",
                "destination": "Đi đường nhánh bên trái về {destination}"
            },
            "slight right": {
                "default": "Đi đường nhánh bên phải",
                "name": "Đi đường nhánh {way_name} bên phải",
                "destination": "Đi đường nhánh bên phải về {destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "Đi vào bùng binh",
                    "name": "Đi vào bùng binh và ra tại {way_name}",
                    "destination": "Đi vào bùng binh và ra về {destination}"
                },
                "name": {
                    "default": "Đi vào {rotary_name}",
                    "name": "Đi vào {rotary_name} và ra tại {way_name}",
                    "destination": "Đi và {rotary_name} và ra về {destination}"
                },
                "exit": {
                    "default": "Đi vào bùng binh và ra tại đường {exit_number}",
                    "name": "Đi vào bùng binh và ra tại đường {exit_number} tức {way_name}",
                    "destination": "Đi vào bùng binh và ra tại đường {exit_number} về {destination}"
                },
                "name_exit": {
                    "default": "Đi vào {rotary_name} và ra tại đường {exit_number}",
                    "name": "Đi vào {rotary_name} và ra tại đường {exit_number} tức {way_name}",
                    "destination": "Đi vào {rotary_name} và ra tại đường {exit_number} về {destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "Đi vào bùng binh và ra tại đường {exit_number}",
                    "name": "Đi vào bùng binh và ra tại đường {exit_number} tức {way_name}",
                    "destination": "Đi vào bùng binh và ra tại đường {exit_number} về {destination}"
                },
                "default": {
                    "default": "Đi vào bùng binh",
                    "name": "Đi vào bùng binh và ra tại {way_name}",
                    "destination": "Đi vào bùng binh và ra về {destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "Quẹo {modifier}",
                "name": "Quẹo {modifier} vào {way_name}",
                "destination": "Quẹo {modifier} về {destination}"
            },
            "left": {
                "default": "Quẹo trái",
                "name": "Quẹo trái vào {way_name}",
                "destination": "Quẹo trái về {destination}"
            },
            "right": {
                "default": "Quẹo phải",
                "name": "Quẹo phải vào {way_name}",
                "destination": "Quẹo phải về {destination}"
            },
            "straight": {
                "default": "Chạy thẳng",
                "name": "Chạy tiếp trên {way_name}",
                "destination": "Chạy tiếp về {destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "Ra bùng binh",
                "name": "Ra bùng binh vào {way_name}",
                "destination": "Ra bùng binh về {destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "Ra bùng binh",
                "name": "Ra bùng binh vào {way_name}",
                "destination": "Ra bùng binh về {destination}"
            }
        },
        "turn": {
            "default": {
                "default": "Quẹo {modifier}",
                "name": "Quẹo {modifier} vào {way_name}",
                "destination": "Quẹo {modifier} về {destination}"
            },
            "left": {
                "default": "Quẹo trái",
                "name": "Quẹo trái vào {way_name}",
                "destination": "Quẹo trái về {destination}"
            },
            "right": {
                "default": "Quẹo phải",
                "name": "Quẹo phải vào {way_name}",
                "destination": "Quẹo phải về {destination}"
            },
            "straight": {
                "default": "Chạy thẳng",
                "name": "Chạy thẳng vào {way_name}",
                "destination": "Chạy thẳng về {destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "Chạy thẳng"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],47:[function(_dereq_,module,exports){
module.exports={
    "meta": {
        "capitalizeFirstLetter": false
    },
    "v5": {
        "constants": {
            "ordinalize": {
                "1": "第一",
                "2": "第二",
                "3": "第三",
                "4": "第四",
                "5": "第五",
                "6": "第六",
                "7": "第七",
                "8": "第八",
                "9": "第九",
                "10": "第十"
            },
            "direction": {
                "north": "北",
                "northeast": "东北",
                "east": "东",
                "southeast": "东南",
                "south": "南",
                "southwest": "西南",
                "west": "西",
                "northwest": "西北"
            },
            "modifier": {
                "left": "向左",
                "right": "向右",
                "sharp left": "急向左",
                "sharp right": "急向右",
                "slight left": "稍向左",
                "slight right": "稍向右",
                "straight": "直行",
                "uturn": "调头"
            },
            "lanes": {
                "xo": "靠右行驶",
                "ox": "靠左行驶",
                "xox": "保持在道路中间行驶",
                "oxo": "保持在道路左侧或右侧行驶"
            }
        },
        "modes": {
            "ferry": {
                "default": "乘坐轮渡",
                "name": "乘坐{way_name}轮渡",
                "destination": "乘坐开往{destination}的轮渡"
            }
        },
        "phrase": {
            "two linked by distance": "{instruction_one}，{distance}后{instruction_two}",
            "two linked": "{instruction_one}，随后{instruction_two}",
            "one in distance": "{distance}后{instruction_one}",
            "name and ref": "{name}（{ref}）",
            "exit with number": "出口{exit}"
        },
        "arrive": {
            "default": {
                "default": "您已经到达您的{nth}个目的地",
                "upcoming": "您即将到达您的{nth}个目的地",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}"
            },
            "left": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路左侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路左侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您左边。"
            },
            "right": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路右侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路右侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您右边。"
            },
            "sharp left": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路左侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路左侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您左边。"
            },
            "sharp right": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路右侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路右侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您右边。"
            },
            "slight right": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路左侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路左侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您右边。"
            },
            "slight left": {
                "default": "您已经到达您的{nth}个目的地，目的地在道路右侧",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在道路右侧",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您左边。"
            },
            "straight": {
                "default": "您已经到达您的{nth}个目的地，目的地在您正前方",
                "upcoming": "您即将到达您的{nth}个目的地，目的地在您正前方",
                "short": "已到达目的地",
                "short-upcoming": "即将到达目的地",
                "named": "您已到达{waypoint_name}，目的地在您前方。"
            }
        },
        "continue": {
            "default": {
                "default": "{modifier}行驶",
                "name": "在{way_name}上继续{modifier}行驶",
                "destination": "{modifier}行驶，{destination}方向",
                "exit": "{modifier}行驶，驶入{way_name}"
            },
            "straight": {
                "default": "继续直行",
                "name": "在{way_name}上继续直行",
                "destination": "继续直行，前往{destination}",
                "distance": "继续直行{distance}",
                "namedistance": "继续在{way_name}上直行{distance}"
            },
            "sharp left": {
                "default": "前方左急转弯",
                "name": "前方左急转弯，继续在{way_name}上行驶",
                "destination": "左急转弯，前往{destination}"
            },
            "sharp right": {
                "default": "前方右急转弯",
                "name": "前方右急转弯，继续在{way_name}上行驶",
                "destination": "右急转弯，前往{destination}"
            },
            "slight left": {
                "default": "前方稍向左转",
                "name": "前方稍向左转，继续在{way_name}上行驶",
                "destination": "稍向左转，前往{destination}"
            },
            "slight right": {
                "default": "前方稍向右转",
                "name": "前方稍向右转，继续在{way_name}上行驶",
                "destination": "前方稍向右转，前往{destination}"
            },
            "uturn": {
                "default": "前方调头",
                "name": "前方调头，继续在{way_name}上行驶",
                "destination": "前方调头，前往{destination}"
            }
        },
        "depart": {
            "default": {
                "default": "出发向{direction}",
                "name": "出发向{direction}，驶入{way_name}",
                "namedistance": "出发向{direction}，在{way_name}上继续行驶{distance}"
            }
        },
        "end of road": {
            "default": {
                "default": "{modifier}行驶",
                "name": "{modifier}行驶，驶入{way_name}",
                "destination": "{modifier}行驶，前往{destination}"
            },
            "straight": {
                "default": "继续直行",
                "name": "继续直行，驶入{way_name}",
                "destination": "继续直行，前往{destination}"
            },
            "uturn": {
                "default": "在道路尽头调头",
                "name": "在道路尽头调头驶入{way_name}",
                "destination": "在道路尽头调头，前往{destination}"
            }
        },
        "fork": {
            "default": {
                "default": "在岔道保持{modifier}",
                "name": "在岔道口保持{modifier}，驶入{way_name}",
                "destination": "在岔道口保持{modifier}，前往{destination}"
            },
            "slight left": {
                "default": "在岔道口保持左侧行驶",
                "name": "在岔道口保持左侧行驶，驶入{way_name}",
                "destination": "在岔道口保持左侧行驶，前往{destination}"
            },
            "slight right": {
                "default": "在岔道口保持右侧行驶",
                "name": "在岔道口保持右侧行驶，驶入{way_name}",
                "destination": "在岔道口保持右侧行驶，前往{destination}"
            },
            "sharp left": {
                "default": "在岔道口左急转弯",
                "name": "在岔道口左急转弯，驶入{way_name}",
                "destination": "在岔道口左急转弯，前往{destination}"
            },
            "sharp right": {
                "default": "在岔道口右急转弯",
                "name": "在岔道口右急转弯，驶入{way_name}",
                "destination": "在岔道口右急转弯，前往{destination}"
            },
            "uturn": {
                "default": "前方调头",
                "name": "前方调头，驶入{way_name}",
                "destination": "前方调头，前往{destination}"
            }
        },
        "merge": {
            "default": {
                "default": "{modifier}并道",
                "name": "{modifier}并道，驶入{way_name}",
                "destination": "{modifier}并道，前往{destination}"
            },
            "straight": {
                "default": "直行并道",
                "name": "直行并道，驶入{way_name}",
                "destination": "直行并道，前往{destination}"
            },
            "slight left": {
                "default": "稍向左并道",
                "name": "稍向左并道，驶入{way_name}",
                "destination": "稍向左并道，前往{destination}"
            },
            "slight right": {
                "default": "稍向右并道",
                "name": "稍向右并道，驶入{way_name}",
                "destination": "稍向右并道，前往{destination}"
            },
            "sharp left": {
                "default": "急向左并道",
                "name": "急向左并道，驶入{way_name}",
                "destination": "急向左并道，前往{destination}"
            },
            "sharp right": {
                "default": "急向右并道",
                "name": "急向右并道，驶入{way_name}",
                "destination": "急向右并道，前往{destination}"
            },
            "uturn": {
                "default": "前方调头",
                "name": "前方调头，驶入{way_name}",
                "destination": "前方调头，前往{destination}"
            }
        },
        "new name": {
            "default": {
                "default": "继续{modifier}",
                "name": "继续{modifier}，驶入{way_name}",
                "destination": "继续{modifier}，前往{destination}"
            },
            "straight": {
                "default": "继续直行",
                "name": "继续在{way_name}上直行",
                "destination": "继续直行，前往{destination}"
            },
            "sharp left": {
                "default": "前方左急转弯",
                "name": "前方左急转弯，驶入{way_name}",
                "destination": "左急转弯，前往{destination}"
            },
            "sharp right": {
                "default": "前方右急转弯",
                "name": "前方右急转弯，驶入{way_name}",
                "destination": "右急转弯，前往{destination}"
            },
            "slight left": {
                "default": "继续稍向左",
                "name": "继续稍向左，驶入{way_name}",
                "destination": "继续稍向左，前往{destination}"
            },
            "slight right": {
                "default": "继续稍向右",
                "name": "继续稍向右，驶入{way_name}",
                "destination": "继续稍向右，前往{destination}"
            },
            "uturn": {
                "default": "前方调头",
                "name": "前方调头，上{way_name}",
                "destination": "前方调头，前往{destination}"
            }
        },
        "notification": {
            "default": {
                "default": "继续{modifier}",
                "name": "继续{modifier}，驶入{way_name}",
                "destination": "继续{modifier}，前往{destination}"
            },
            "uturn": {
                "default": "前方调头",
                "name": "前方调头，驶入{way_name}",
                "destination": "前方调头，前往{destination}"
            }
        },
        "off ramp": {
            "default": {
                "default": "下匝道",
                "name": "下匝道，驶入{way_name}",
                "destination": "下匝道，前往{destination}",
                "exit": "从{exit}出口驶出",
                "exit_destination": "从{exit}出口驶出，前往{destination}"
            },
            "left": {
                "default": "下左侧匝道",
                "name": "下左侧匝道，上{way_name}",
                "destination": "下左侧匝道，前往{destination}",
                "exit": "从左侧{exit}出口驶出",
                "exit_destination": "从左侧{exit}出口驶出，前往{destination}"
            },
            "right": {
                "default": "下右侧匝道",
                "name": "下右侧匝道，驶入{way_name}",
                "destination": "下右侧匝道，前往{destination}",
                "exit": "从右侧{exit}出口驶出",
                "exit_destination": "从右侧{exit}出口驶出，前往{destination}"
            },
            "sharp left": {
                "default": "急向左下匝道",
                "name": "急向左下匝道，驶入{way_name}",
                "destination": "急向左下匝道，前往{destination}",
                "exit": "从左侧{exit}出口驶出",
                "exit_destination": "从左侧{exit}出口驶出，前往{destination}"
            },
            "sharp right": {
                "default": "急向右下匝道",
                "name": "急向右下匝道，驶入{way_name}",
                "destination": "急向右下匝道，前往{destination}",
                "exit": "从右侧{exit}出口驶出",
                "exit_destination": "从右侧{exit}出口驶出，前往{destination}"
            },
            "slight left": {
                "default": "稍向左下匝道",
                "name": "稍向左下匝道，驶入{way_name}",
                "destination": "稍向左下匝道，前往{destination}",
                "exit": "从左侧{exit}出口驶出",
                "exit_destination": "从左侧{exit}出口驶出，前往{destination}"
            },
            "slight right": {
                "default": "稍向右下匝道",
                "name": "稍向右下匝道，驶入{way_name}",
                "destination": "稍向右下匝道，前往{destination}",
                "exit": "从右侧{exit}出口驶出",
                "exit_destination": "从右侧{exit}出口驶出，前往{destination}"
            }
        },
        "on ramp": {
            "default": {
                "default": "上匝道",
                "name": "上匝道，驶入{way_name}",
                "destination": "上匝道，前往{destination}"
            },
            "left": {
                "default": "上左侧匝道",
                "name": "上左侧匝道，驶入{way_name}",
                "destination": "上左侧匝道，前往{destination}"
            },
            "right": {
                "default": "上右侧匝道",
                "name": "上右侧匝道，驶入{way_name}",
                "destination": "上右侧匝道，前往{destination}"
            },
            "sharp left": {
                "default": "急向左上匝道",
                "name": "急向左上匝道，驶入{way_name}",
                "destination": "急向左上匝道，前往{destination}"
            },
            "sharp right": {
                "default": "急向右上匝道",
                "name": "急向右上匝道，驶入{way_name}",
                "destination": "急向右上匝道，前往{destination}"
            },
            "slight left": {
                "default": "稍向左上匝道",
                "name": "稍向左上匝道，驶入{way_name}",
                "destination": "稍向左上匝道，前往{destination}"
            },
            "slight right": {
                "default": "稍向右上匝道",
                "name": "稍向右上匝道，驶入{way_name}",
                "destination": "稍向右上匝道，前往{destination}"
            }
        },
        "rotary": {
            "default": {
                "default": {
                    "default": "进入环岛",
                    "name": "通过环岛后驶入{way_name}",
                    "destination": "通过环岛后前往{destination}"
                },
                "name": {
                    "default": "进入{rotary_name}环岛",
                    "name": "通过{rotary_name}环岛后驶入{way_name}",
                    "destination": "通过{rotary_name}环岛后前往{destination}"
                },
                "exit": {
                    "default": "进入环岛后从{exit_number}出口驶出",
                    "name": "进入环岛后从{exit_number}出口驶出，上{way_name}",
                    "destination": "进入环岛后从{exit_number}出口驶出，前往{destination}"
                },
                "name_exit": {
                    "default": "进入{rotary_name}环岛后从{exit_number}出口驶出",
                    "name": "进入{rotary_name}环岛后从{exit_number}出口驶出，上{way_name}",
                    "destination": "进入{rotary_name}环岛后从{exit_number}出口驶出，前往{destination}"
                }
            }
        },
        "roundabout": {
            "default": {
                "exit": {
                    "default": "进入环岛后从{exit_number}出口驶出",
                    "name": "进入环岛后从{exit_number}出口驶出，上{way_name}",
                    "destination": "进入环岛后从{exit_number}出口驶出，前往{destination}"
                },
                "default": {
                    "default": "进入环岛",
                    "name": "通过环岛后驶入{way_name}",
                    "destination": "通过环岛后前往{destination}"
                }
            }
        },
        "roundabout turn": {
            "default": {
                "default": "{modifier}转弯",
                "name": "{modifier}转弯，驶入{way_name}",
                "destination": "{modifier}转弯，前往{destination}"
            },
            "left": {
                "default": "左转",
                "name": "左转，驶入{way_name}",
                "destination": "左转，前往{destination}"
            },
            "right": {
                "default": "右转",
                "name": "右转，驶入{way_name}",
                "destination": "右转，前往{destination}"
            },
            "straight": {
                "default": "继续直行",
                "name": "继续直行，驶入{way_name}",
                "destination": "继续直行，前往{destination}"
            }
        },
        "exit roundabout": {
            "default": {
                "default": "驶离环岛",
                "name": "驶离环岛，驶入{way_name}",
                "destination": "驶离环岛，前往{destination}"
            }
        },
        "exit rotary": {
            "default": {
                "default": "驶离环岛",
                "name": "驶离环岛，驶入{way_name}",
                "destination": "驶离环岛，前往{destination}"
            }
        },
        "turn": {
            "default": {
                "default": "{modifier}转弯",
                "name": "{modifier}转弯，驶入{way_name}",
                "destination": "{modifier}转弯，前往{destination}"
            },
            "left": {
                "default": "左转",
                "name": "左转，驶入{way_name}",
                "destination": "左转，前往{destination}"
            },
            "right": {
                "default": "右转",
                "name": "右转，驶入{way_name}",
                "destination": "右转，前往{destination}"
            },
            "straight": {
                "default": "直行",
                "name": "直行，驶入{way_name}",
                "destination": "直行，前往{destination}"
            }
        },
        "use lane": {
            "no_lanes": {
                "default": "继续直行"
            },
            "default": {
                "default": "{lane_instruction}"
            }
        }
    }
}

},{}],48:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	module.exports = L.Class.extend({
		options: {
			timeout: 500,
			blurTimeout: 100,
			noResultsMessage: 'No results found.'
		},

		initialize: function(elem, callback, context, options) {
			L.setOptions(this, options);

			this._elem = elem;
			this._resultFn = options.resultFn ? L.Util.bind(options.resultFn, options.resultContext) : null;
			this._autocomplete = options.autocompleteFn ? L.Util.bind(options.autocompleteFn, options.autocompleteContext) : null;
			this._selectFn = L.Util.bind(callback, context);
			this._container = L.DomUtil.create('div', 'leaflet-routing-geocoder-result');
			this._resultTable = L.DomUtil.create('table', '', this._container);

			// TODO: looks a bit like a kludge to register same for input and keypress -
			// browsers supporting both will get duplicate events; just registering
			// input will not catch enter, though.
			L.DomEvent.addListener(this._elem, 'input', this._keyPressed, this);
			L.DomEvent.addListener(this._elem, 'keypress', this._keyPressed, this);
			L.DomEvent.addListener(this._elem, 'keydown', this._keyDown, this);
			L.DomEvent.addListener(this._elem, 'blur', function() {
				if (this._isOpen) {
					this.close();
				}
			}, this);
		},

		close: function() {
			L.DomUtil.removeClass(this._container, 'leaflet-routing-geocoder-result-open');
			this._isOpen = false;
		},

		_open: function() {
			var rect = this._elem.getBoundingClientRect();
			if (!this._container.parentElement) {
				// See notes section under https://developer.mozilla.org/en-US/docs/Web/API/Window/scrollX
				// This abomination is required to support all flavors of IE
				var scrollX = (window.pageXOffset !== undefined) ? window.pageXOffset
					: (document.documentElement || document.body.parentNode || document.body).scrollLeft;
				var scrollY = (window.pageYOffset !== undefined) ? window.pageYOffset
					: (document.documentElement || document.body.parentNode || document.body).scrollTop;
				this._container.style.left = (rect.left + scrollX) + 'px';
				this._container.style.top = (rect.bottom + scrollY) + 'px';
				this._container.style.width = (rect.right - rect.left) + 'px';
				document.body.appendChild(this._container);
			}

			L.DomUtil.addClass(this._container, 'leaflet-routing-geocoder-result-open');
			this._isOpen = true;
		},

		_setResults: function(results) {
			var i,
			    tr,
			    td,
			    text;

			delete this._selection;
			this._results = results;

			while (this._resultTable.firstChild) {
				this._resultTable.removeChild(this._resultTable.firstChild);
			}

			for (i = 0; i < results.length; i++) {
				tr = L.DomUtil.create('tr', '', this._resultTable);
				tr.setAttribute('data-result-index', i);
				td = L.DomUtil.create('td', '', tr);
				text = document.createTextNode(results[i].name);
				td.appendChild(text);
				// mousedown + click because:
				// http://stackoverflow.com/questions/10652852/jquery-fire-click-before-blur-event
				L.DomEvent.addListener(td, 'mousedown', L.DomEvent.preventDefault);
				L.DomEvent.addListener(td, 'click', this._createClickListener(results[i]));
			}

			if (!i) {
				tr = L.DomUtil.create('tr', '', this._resultTable);
				td = L.DomUtil.create('td', 'leaflet-routing-geocoder-no-results', tr);
				td.innerHTML = this.options.noResultsMessage;
			}

			this._open();

			if (results.length > 0) {
				// Select the first entry
				this._select(1);
			}
		},

		_createClickListener: function(r) {
			var resultSelected = this._resultSelected(r);
			return L.bind(function() {
				this._elem.blur();
				resultSelected();
			}, this);
		},

		_resultSelected: function(r) {
			return L.bind(function() {
				this.close();
				this._elem.value = r.name;
				this._lastCompletedText = r.name;
				this._selectFn(r);
			}, this);
		},

		_keyPressed: function(e) {
			var index;

			if (this._isOpen && e.keyCode === 13 && this._selection) {
				index = parseInt(this._selection.getAttribute('data-result-index'), 10);
				this._resultSelected(this._results[index])();
				L.DomEvent.preventDefault(e);
				return;
			}

			if (e.keyCode === 13) {
				L.DomEvent.preventDefault(e);
				this._complete(this._resultFn, true);
				return;
			}

			if (this._autocomplete && document.activeElement === this._elem) {
				if (this._timer) {
					clearTimeout(this._timer);
				}
				this._timer = setTimeout(L.Util.bind(function() { this._complete(this._autocomplete); }, this),
					this.options.timeout);
				return;
			}

			this._unselect();
		},

		_select: function(dir) {
			var sel = this._selection;
			if (sel) {
				L.DomUtil.removeClass(sel.firstChild, 'leaflet-routing-geocoder-selected');
				sel = sel[dir > 0 ? 'nextSibling' : 'previousSibling'];
			}
			if (!sel) {
				sel = this._resultTable[dir > 0 ? 'firstChild' : 'lastChild'];
			}

			if (sel) {
				L.DomUtil.addClass(sel.firstChild, 'leaflet-routing-geocoder-selected');
				this._selection = sel;
			}
		},

		_unselect: function() {
			if (this._selection) {
				L.DomUtil.removeClass(this._selection.firstChild, 'leaflet-routing-geocoder-selected');
			}
			delete this._selection;
		},

		_keyDown: function(e) {
			if (this._isOpen) {
				switch (e.keyCode) {
				// Escape
				case 27:
					this.close();
					L.DomEvent.preventDefault(e);
					return;
				// Up
				case 38:
					this._select(-1);
					L.DomEvent.preventDefault(e);
					return;
				// Down
				case 40:
					this._select(1);
					L.DomEvent.preventDefault(e);
					return;
				}
			}
		},

		_complete: function(completeFn, trySelect) {
			var v = this._elem.value;
			function completeResults(results) {
				this._lastCompletedText = v;
				if (trySelect && results.length === 1) {
					this._resultSelected(results[0])();
				} else {
					this._setResults(results);
				}
			}

			if (!v) {
				return;
			}

			if (v !== this._lastCompletedText) {
				completeFn(v, completeResults, this);
			} else if (trySelect) {
				completeResults.call(this, this._results);
			}
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],49:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	var Itinerary = _dereq_('./itinerary');
	var Line = _dereq_('./line');
	var Plan = _dereq_('./plan');
	var OSRMv1 = _dereq_('./osrm-v1');

	module.exports = Itinerary.extend({
		options: {
			fitSelectedRoutes: 'smart',
			routeLine: function(route, options) { return new Line(route, options); },
			autoRoute: true,
			routeWhileDragging: false,
			routeDragInterval: 500,
			waypointMode: 'connect',
			showAlternatives: false,
			defaultErrorHandler: function(e) {
				console.error('Routing error:', e.error);
			}
		},

		initialize: function(options) {
			L.Util.setOptions(this, options);

			this._router = this.options.router || new OSRMv1(options);
			this._plan = this.options.plan || new Plan(this.options.waypoints, options);
			this._requestCount = 0;

			Itinerary.prototype.initialize.call(this, options);

			this.on('routeselected', this._routeSelected, this);
			if (this.options.defaultErrorHandler) {
				this.on('routingerror', this.options.defaultErrorHandler);
			}
			this._plan.on('waypointschanged', this._onWaypointsChanged, this);
			if (options.routeWhileDragging) {
				this._setupRouteDragging();
			}
				if (this.options.autoRoute) {
				this.route();
			}
		},

		_onZoomEnd: function() {
			if (!this._selectedRoute ||
				!this._router.requiresMoreDetail) {
				return;
			}

			var map = this._map;
			if (this._router.requiresMoreDetail(this._selectedRoute,
					map.getZoom(), map.getBounds())) {
				this.route({
					callback: L.bind(function(err, routes) {
						var i;
						if (!err) {
							for (i = 0; i < routes.length; i++) {
								this._routes[i].properties = routes[i].properties;
							}
							this._updateLineCallback(err, routes);
						}

					}, this),
					simplifyGeometry: false,
					geometryOnly: true
				});
			}
		},

		onAdd: function(map) {
			if (this.options.autoRoute) {
				this.route();
			}

			var container = Itinerary.prototype.onAdd.call(this, map);

			this._map = map;
			this._map.addLayer(this._plan);

			this._map.on('zoomend', this._onZoomEnd, this);

			if (this._plan.options.geocoder) {
				container.insertBefore(this._plan.createGeocoders(), container.firstChild);
			}

			return container;
		},

		onRemove: function(map) {
			map.off('zoomend', this._onZoomEnd, this);
			if (this._line) {
				map.removeLayer(this._line);
			}
			map.removeLayer(this._plan);
			if (this._alternatives && this._alternatives.length > 0) {
				for (var i = 0, len = this._alternatives.length; i < len; i++) {
					map.removeLayer(this._alternatives[i]);
				}
			}
			return Itinerary.prototype.onRemove.call(this, map);
		},

		getWaypoints: function() {
			return this._plan.getWaypoints();
		},

		setWaypoints: function(waypoints) {
			this._plan.setWaypoints(waypoints);
			return this;
		},

		spliceWaypoints: function() {
			var removed = this._plan.spliceWaypoints.apply(this._plan, arguments);
			return removed;
		},

		getPlan: function() {
			return this._plan;
		},

		getRouter: function() {
			return this._router;
		},

		_routeSelected: function(e) {
			var route = this._selectedRoute = e.route,
				alternatives = this.options.showAlternatives && e.alternatives,
				fitMode = this.options.fitSelectedRoutes,
				fitBounds =
					(fitMode === 'smart' && !this._waypointsVisible()) ||
					(fitMode !== 'smart' && fitMode);

			this._updateLines({route: route, alternatives: alternatives});

			if (fitBounds) {
				this._map.fitBounds(this._line.getBounds());
			}

			if (this.options.waypointMode === 'snap') {
				this._plan.off('waypointschanged', this._onWaypointsChanged, this);
				this.setWaypoints(route.waypoints);
				this._plan.on('waypointschanged', this._onWaypointsChanged, this);
			}
		},

		_waypointsVisible: function() {
			var wps = this.getWaypoints(),
				mapSize,
				bounds,
				boundsSize,
				i,
				p;

			try {
				mapSize = this._map.getSize();

				for (i = 0; i < wps.length; i++) {
					p = this._map.latLngToLayerPoint(wps[i].latLng);

					if (bounds) {
						bounds.extend(p);
					} else {
						bounds = L.bounds([p]);
					}
				}

				boundsSize = bounds.getSize();
				return (boundsSize.x > mapSize.x / 5 ||
					boundsSize.y > mapSize.y / 5) && this._waypointsInViewport();

			} catch (e) {
				return false;
			}
		},

		_waypointsInViewport: function() {
			var wps = this.getWaypoints(),
				mapBounds,
				i;

			try {
				mapBounds = this._map.getBounds();
			} catch (e) {
				return false;
			}

			for (i = 0; i < wps.length; i++) {
				if (mapBounds.contains(wps[i].latLng)) {
					return true;
				}
			}

			return false;
		},

		_updateLines: function(routes) {
			var addWaypoints = this.options.addWaypoints !== undefined ?
				this.options.addWaypoints : true;
			this._clearLines();

			// add alternatives first so they lie below the main route
			this._alternatives = [];
			if (routes.alternatives) routes.alternatives.forEach(function(alt, i) {
				this._alternatives[i] = this.options.routeLine(alt,
					L.extend({
						isAlternative: true
					}, this.options.altLineOptions || this.options.lineOptions));
				this._alternatives[i].addTo(this._map);
				this._hookAltEvents(this._alternatives[i]);
			}, this);

			this._line = this.options.routeLine(routes.route,
				L.extend({
					addWaypoints: addWaypoints,
					extendToWaypoints: this.options.waypointMode === 'connect'
				}, this.options.lineOptions));
			this._line.addTo(this._map);
			this._hookEvents(this._line);
		},

		_hookEvents: function(l) {
			l.on('linetouched', function(e) {
				this._plan.dragNewWaypoint(e);
			}, this);
		},

		_hookAltEvents: function(l) {
			l.on('linetouched', function(e) {
				var alts = this._routes.slice();
				var selected = alts.splice(e.target._route.routesIndex, 1)[0];
				this.fire('routeselected', {route: selected, alternatives: alts});
			}, this);
		},

		_onWaypointsChanged: function(e) {
			if (this.options.autoRoute) {
				this.route({});
			}
			if (!this._plan.isReady()) {
				this._clearLines();
				this._clearAlts();
			}
			this.fire('waypointschanged', {waypoints: e.waypoints});
		},

		_setupRouteDragging: function() {
			var timer = 0,
				waypoints;

			this._plan.on('waypointdrag', L.bind(function(e) {
				waypoints = e.waypoints;

				if (!timer) {
					timer = setTimeout(L.bind(function() {
						this.route({
							waypoints: waypoints,
							geometryOnly: true,
							callback: L.bind(this._updateLineCallback, this)
						});
						timer = undefined;
					}, this), this.options.routeDragInterval);
				}
			}, this));
			this._plan.on('waypointdragend', function() {
				if (timer) {
					clearTimeout(timer);
					timer = undefined;
				}
				this.route();
			}, this);
		},

		_updateLineCallback: function(err, routes) {
			if (!err) {
				routes = routes.slice();
				var selected = routes.splice(this._selectedRoute.routesIndex, 1)[0];
				this._updateLines({
					route: selected,
					alternatives: this.options.showAlternatives ? routes : []
				});
			} else if (err.type !== 'abort') {
				this._clearLines();
			}
		},

		route: function(options) {
			var ts = ++this._requestCount,
				wps;

			if (this._pendingRequest && this._pendingRequest.abort) {
				this._pendingRequest.abort();
				this._pendingRequest = null;
			}

			options = options || {};

			if (this._plan.isReady()) {
				if (this.options.useZoomParameter) {
					options.z = this._map && this._map.getZoom();
				}

				wps = options && options.waypoints || this._plan.getWaypoints();
				this.fire('routingstart', {waypoints: wps});
				this._pendingRequest = this._router.route(wps, function(err, routes) {
					this._pendingRequest = null;

					if (options.callback) {
						return options.callback.call(this, err, routes);
					}

					// Prevent race among multiple requests,
					// by checking the current request's count
					// against the last request's; ignore result if
					// this isn't the last request.
					if (ts === this._requestCount) {
						this._clearLines();
						this._clearAlts();
						if (err && err.type !== 'abort') {
							this.fire('routingerror', {error: err});
							return;
						}

						routes.forEach(function(route, i) { route.routesIndex = i; });

						if (!options.geometryOnly) {
							this.fire('routesfound', {waypoints: wps, routes: routes});
							this.setAlternatives(routes);
						} else {
							var selectedRoute = routes.splice(0,1)[0];
							this._routeSelected({route: selectedRoute, alternatives: routes});
						}
					}
				}, this, options);
			}
		},

		_clearLines: function() {
			if (this._line) {
				this._map.removeLayer(this._line);
				delete this._line;
			}
			if (this._alternatives && this._alternatives.length) {
				for (var i in this._alternatives) {
					this._map.removeLayer(this._alternatives[i]);
				}
				this._alternatives = [];
			}
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./itinerary":55,"./line":56,"./osrm-v1":59,"./plan":60}],50:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	module.exports = L.Control.extend({
		options: {
			header: 'Routing error',
			formatMessage: function(error) {
				if (error.status < 0) {
					return 'Calculating the route caused an error. Technical description follows: <code><pre>' +
						error.message + '</pre></code';
				} else {
					return 'The route could not be calculated. ' +
						error.message;
				}
			}
		},

		initialize: function(routingControl, options) {
			L.Control.prototype.initialize.call(this, options);
			routingControl
				.on('routingerror', L.bind(function(e) {
					if (this._element) {
						this._element.children[1].innerHTML = this.options.formatMessage(e.error);
						this._element.style.visibility = 'visible';
					}
				}, this))
				.on('routingstart', L.bind(function() {
					if (this._element) {
						this._element.style.visibility = 'hidden';
					}
				}, this));
		},

		onAdd: function() {
			var header,
				message;

			this._element = L.DomUtil.create('div', 'leaflet-bar leaflet-routing-error');
			this._element.style.visibility = 'hidden';

			header = L.DomUtil.create('h3', null, this._element);
			message = L.DomUtil.create('span', null, this._element);

			header.innerHTML = this.options.header;

			return this._element;
		},

		onRemove: function() {
			delete this._element;
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],51:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	var Localization = _dereq_('./localization');

	module.exports = L.Class.extend({
		options: {
			units: 'metric',
			unitNames: null,
			language: 'en',
			roundingSensitivity: 1,
			distanceTemplate: '{value} {unit}'
		},

		initialize: function(options) {
			L.setOptions(this, options);

			var langs = L.Util.isArray(this.options.language) ?
				this.options.language :
				[this.options.language, 'en'];
			this._localization = new Localization(langs);
		},

		formatDistance: function(d /* Number (meters) */, sensitivity) {
			var un = this.options.unitNames || this._localization.localize('units'),
				simpleRounding = sensitivity <= 0,
				round = simpleRounding ? function(v) { return v; } : L.bind(this._round, this),
			    v,
			    yards,
				data,
				pow10;

			if (this.options.units === 'imperial') {
				yards = d / 0.9144;
				if (yards >= 1000) {
					data = {
						value: round(d / 1609.344, sensitivity),
						unit: un.miles
					};
				} else {
					data = {
						value: round(yards, sensitivity),
						unit: un.yards
					};
				}
			} else {
				v = round(d, sensitivity);
				data = {
					value: v >= 1000 ? (v / 1000) : v,
					unit: v >= 1000 ? un.kilometers : un.meters
				};
			}

			if (simpleRounding) {
				data.value = data.value.toFixed(-sensitivity);
			}

			return L.Util.template(this.options.distanceTemplate, data);
		},

		_round: function(d, sensitivity) {
			var s = sensitivity || this.options.roundingSensitivity,
				pow10 = Math.pow(10, (Math.floor(d / s) + '').length - 1),
				r = Math.floor(d / pow10),
				p = (r > 5) ? pow10 : pow10 / 2;

			return Math.round(d / p) * p;
		},

		formatTime: function(t /* Number (seconds) */) {
			var un = this.options.unitNames || this._localization.localize('units');
			// More than 30 seconds precision looks ridiculous
			t = Math.round(t / 30) * 30;

			if (t > 86400) {
				return Math.round(t / 3600) + ' ' + un.hours;
			} else if (t > 3600) {
				return Math.floor(t / 3600) + ' ' + un.hours + ' ' +
					Math.round((t % 3600) / 60) + ' ' + un.minutes;
			} else if (t > 300) {
				return Math.round(t / 60) + ' ' + un.minutes;
			} else if (t > 60) {
				return Math.floor(t / 60) + ' ' + un.minutes +
					(t % 60 !== 0 ? ' ' + (t % 60) + ' ' + un.seconds : '');
			} else {
				return t + ' ' + un.seconds;
			}
		},

		formatInstruction: function(instr, i) {
			if (instr.text === undefined) {
				return this.capitalize(L.Util.template(this._getInstructionTemplate(instr, i),
					L.extend({}, instr, {
						exitStr: instr.exit ? this._localization.localize('formatOrder')(instr.exit) : '',
						dir: this._localization.localize(['directions', instr.direction]),
						modifier: this._localization.localize(['directions', instr.modifier])
					})));
			} else {
				return instr.text;
			}
		},

		getIconName: function(instr, i) {
			switch (instr.type) {
			case 'Head':
				if (i === 0) {
					return 'depart';
				}
				break;
			case 'WaypointReached':
				return 'via';
			case 'Roundabout':
				return 'enter-roundabout';
			case 'DestinationReached':
				return 'arrive';
			}

			switch (instr.modifier) {
			case 'Straight':
				return 'continue';
			case 'SlightRight':
				return 'bear-right';
			case 'Right':
				return 'turn-right';
			case 'SharpRight':
				return 'sharp-right';
			case 'TurnAround':
			case 'Uturn':
				return 'u-turn';
			case 'SharpLeft':
				return 'sharp-left';
			case 'Left':
				return 'turn-left';
			case 'SlightLeft':
				return 'bear-left';
			}
		},

		capitalize: function(s) {
			return s.charAt(0).toUpperCase() + s.substring(1);
		},

		_getInstructionTemplate: function(instr, i) {
			var type = instr.type === 'Straight' ? (i === 0 ? 'Head' : 'Continue') : instr.type,
				strings = this._localization.localize(['instructions', type]);

			if (!strings) {
				strings = [
					this._localization.localize(['directions', type]),
					' ' + this._localization.localize(['instructions', 'Onto'])
				];
			}

			return strings[0] + (strings.length > 1 && instr.road ? strings[1] : '');
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./localization":57}],52:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);
	var Autocomplete = _dereq_('./autocomplete');
	var Localization = _dereq_('./localization');

	function selectInputText(input) {
		if (input.setSelectionRange) {
			// On iOS, select() doesn't work
			input.setSelectionRange(0, 9999);
		} else {
			// On at least IE8, setSeleectionRange doesn't exist
			input.select();
		}
	}

	module.exports = L.Class.extend({
		includes: ((typeof L.Evented !== 'undefined' && L.Evented.prototype) || L.Mixin.Events),

		options: {
			createGeocoder: function(i, nWps, options) {
				var container = L.DomUtil.create('div', 'leaflet-routing-geocoder'),
					input = L.DomUtil.create('input', '', container),
					remove = options.addWaypoints ? L.DomUtil.create('span', 'leaflet-routing-remove-waypoint', container) : undefined;

				input.disabled = !options.addWaypoints;

				return {
					container: container,
					input: input,
					closeButton: remove
				};
			},
			geocoderPlaceholder: function(i, numberWaypoints, geocoderElement) {
				var l = new Localization(geocoderElement.options.language).localize('ui');
				return i === 0 ?
					l.startPlaceholder :
					(i < numberWaypoints - 1 ?
						L.Util.template(l.viaPlaceholder, {viaNumber: i}) :
						l.endPlaceholder);
			},

			geocoderClass: function() {
				return '';
			},

			waypointNameFallback: function(latLng) {
				var ns = latLng.lat < 0 ? 'S' : 'N',
					ew = latLng.lng < 0 ? 'W' : 'E',
					lat = (Math.round(Math.abs(latLng.lat) * 10000) / 10000).toString(),
					lng = (Math.round(Math.abs(latLng.lng) * 10000) / 10000).toString();
				return ns + lat + ', ' + ew + lng;
			},
			maxGeocoderTolerance: 200,
			autocompleteOptions: {},
			language: 'en',
		},

		initialize: function(wp, i, nWps, options) {
			L.setOptions(this, options);

			var g = this.options.createGeocoder(i, nWps, this.options),
				closeButton = g.closeButton,
				geocoderInput = g.input;
			geocoderInput.setAttribute('placeholder', this.options.geocoderPlaceholder(i, nWps, this));
			geocoderInput.className = this.options.geocoderClass(i, nWps);

			this._element = g;
			this._waypoint = wp;

			this.update();
			// This has to be here, or geocoder's value will not be properly
			// initialized.
			// TODO: look into why and make _updateWaypointName fix this.
			geocoderInput.value = wp.name;

			L.DomEvent.addListener(geocoderInput, 'click', function() {
				selectInputText(this);
			}, geocoderInput);

			if (closeButton) {
				L.DomEvent.addListener(closeButton, 'click', function() {
					this.fire('delete', { waypoint: this._waypoint });
				}, this);
			}

			new Autocomplete(geocoderInput, function(r) {
					geocoderInput.value = r.name;
					wp.name = r.name;
					wp.latLng = r.center;
					this.fire('geocoded', { waypoint: wp, value: r });
				}, this, L.extend({
					resultFn: this.options.geocoder.geocode,
					resultContext: this.options.geocoder,
					autocompleteFn: this.options.geocoder.suggest,
					autocompleteContext: this.options.geocoder
				}, this.options.autocompleteOptions));
		},

		getContainer: function() {
			return this._element.container;
		},

		setValue: function(v) {
			this._element.input.value = v;
		},

		update: function(force) {
			var wp = this._waypoint,
				wpCoords;

			wp.name = wp.name || '';

			if (wp.latLng && (force || !wp.name)) {
				wpCoords = this.options.waypointNameFallback(wp.latLng);
				if (this.options.geocoder && this.options.geocoder.reverse) {
					this.options.geocoder.reverse(wp.latLng, 67108864 /* zoom 18 */, function(rs) {
						if (rs.length > 0 && rs[0].center.distanceTo(wp.latLng) < this.options.maxGeocoderTolerance) {
							wp.name = rs[0].name;
						} else {
							wp.name = wpCoords;
						}
						this._update();
					}, this);
				} else {
					wp.name = wpCoords;
					this._update();
				}
			}
		},

		focus: function() {
			var input = this._element.input;
			input.focus();
			selectInputText(input);
		},

		_update: function() {
			var wp = this._waypoint,
			    value = wp && wp.name ? wp.name : '';
			this.setValue(value);
			this.fire('reversegeocoded', {waypoint: wp, value: value});
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./autocomplete":48,"./localization":57}],53:[function(_dereq_,module,exports){
(function (global){
var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null),
    Control = _dereq_('./control'),
    Itinerary = _dereq_('./itinerary'),
    Line = _dereq_('./line'),
    OSRMv1 = _dereq_('./osrm-v1'),
    Plan = _dereq_('./plan'),
    Waypoint = _dereq_('./waypoint'),
    Autocomplete = _dereq_('./autocomplete'),
    Formatter = _dereq_('./formatter'),
    GeocoderElement = _dereq_('./geocoder-element'),
    Localization = _dereq_('./localization'),
    ItineraryBuilder = _dereq_('./itinerary-builder'),
    Mapbox = _dereq_('./mapbox'),
    ErrorControl = _dereq_('./error-control');

L.routing = {
    control: function(options) { return new Control(options); },
    itinerary: function(options) {
        return Itinerary(options);
    },
    line: function(route, options) {
        return new Line(route, options);
    },
    plan: function(waypoints, options) {
        return new Plan(waypoints, options);
    },
    waypoint: function(latLng, name, options) {
        return new Waypoint(latLng, name, options);
    },
    osrmv1: function(options) {
        return new OSRMv1(options);
    },
    localization: function(options) {
        return new Localization(options);
    },
    formatter: function(options) {
        return new Formatter(options);
    },
    geocoderElement: function(wp, i, nWps, plan) {
        return new L.Routing.GeocoderElement(wp, i, nWps, plan);
    },
    itineraryBuilder: function(options) {
        return new ItineraryBuilder(options);
    },
    mapbox: function(accessToken, options) {
        return new Mapbox(accessToken, options);
    },
    errorControl: function(routingControl, options) {
        return new ErrorControl(routingControl, options);
    },
    autocomplete: function(elem, callback, context, options) {
        return new Autocomplete(elem, callback, context, options);
    }
};

module.exports = L.Routing = {
    Control: Control,
    Itinerary: Itinerary,
    Line: Line,
    OSRMv1: OSRMv1,
    Plan: Plan,
    Waypoint: Waypoint,
    Autocomplete: Autocomplete,
    Formatter: Formatter,
    GeocoderElement: GeocoderElement,
    Localization: Localization,
    Formatter: Formatter,
    ItineraryBuilder: ItineraryBuilder,

    // Legacy; remove these in next major release
    control: L.routing.control,
    itinerary: L.routing.itinerary,
    line: L.routing.line,
    plan: L.routing.plan,
    waypoint: L.routing.waypoint,
    osrmv1: L.routing.osrmv1,
    geocoderElement: L.routing.geocoderElement,
    mapbox: L.routing.mapbox,
    errorControl: L.routing.errorControl,
};

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./autocomplete":48,"./control":49,"./error-control":50,"./formatter":51,"./geocoder-element":52,"./itinerary":55,"./itinerary-builder":54,"./line":56,"./localization":57,"./mapbox":58,"./osrm-v1":59,"./plan":60,"./waypoint":61}],54:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	module.exports = L.Class.extend({
		options: {
			containerClassName: ''
		},

		initialize: function(options) {
			L.setOptions(this, options);
		},

		createContainer: function(className) {
			var table = L.DomUtil.create('table', (className || '') + ' ' + this.options.containerClassName),
				colgroup = L.DomUtil.create('colgroup', '', table);

			L.DomUtil.create('col', 'leaflet-routing-instruction-icon', colgroup);
			L.DomUtil.create('col', 'leaflet-routing-instruction-text', colgroup);
			L.DomUtil.create('col', 'leaflet-routing-instruction-distance', colgroup);

			return table;
		},

		createStepsContainer: function() {
			return L.DomUtil.create('tbody', '');
		},

		createStep: function(text, distance, icon, steps) {
			var row = L.DomUtil.create('tr', '', steps),
				span,
				td;
			td = L.DomUtil.create('td', '', row);
			span = L.DomUtil.create('span', 'leaflet-routing-icon leaflet-routing-icon-'+icon, td);
			td.appendChild(span);
			td = L.DomUtil.create('td', '', row);
			td.appendChild(document.createTextNode(text));
			td = L.DomUtil.create('td', '', row);
			td.appendChild(document.createTextNode(distance));
			return row;
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],55:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);
	var Formatter = _dereq_('./formatter');
	var ItineraryBuilder = _dereq_('./itinerary-builder');

	module.exports = L.Control.extend({
		includes: ((typeof L.Evented !== 'undefined' && L.Evented.prototype) || L.Mixin.Events),

		options: {
			pointMarkerStyle: {
				radius: 5,
				color: '#03f',
				fillColor: 'white',
				opacity: 1,
				fillOpacity: 0.7
			},
			summaryTemplate: '<h2>{name}</h2><h3>{distance}, {time}</h3>',
			timeTemplate: '{time}',
			containerClassName: '',
			alternativeClassName: '',
			minimizedClassName: '',
			itineraryClassName: '',
			totalDistanceRoundingSensitivity: -1,
			show: true,
			collapsible: undefined,
			collapseBtn: function(itinerary) {
				var collapseBtn = L.DomUtil.create('span', itinerary.options.collapseBtnClass);
				L.DomEvent.on(collapseBtn, 'click', itinerary._toggle, itinerary);
				itinerary._container.insertBefore(collapseBtn, itinerary._container.firstChild);
			},
			collapseBtnClass: 'leaflet-routing-collapse-btn'
		},

		initialize: function(options) {
			L.setOptions(this, options);
			this._formatter = this.options.formatter || new Formatter(this.options);
			this._itineraryBuilder = this.options.itineraryBuilder || new ItineraryBuilder({
				containerClassName: this.options.itineraryClassName
			});
		},

		onAdd: function(map) {
			var collapsible = this.options.collapsible;

			collapsible = collapsible || (collapsible === undefined && map.getSize().x <= 640);

			this._container = L.DomUtil.create('div', 'leaflet-routing-container leaflet-bar ' +
				(!this.options.show ? 'leaflet-routing-container-hide ' : '') +
				(collapsible ? 'leaflet-routing-collapsible ' : '') +
				this.options.containerClassName);
			this._altContainer = this.createAlternativesContainer();
			this._container.appendChild(this._altContainer);
			L.DomEvent.disableClickPropagation(this._container);
			L.DomEvent.addListener(this._container, 'mousewheel', function(e) {
				L.DomEvent.stopPropagation(e);
			});

			if (collapsible) {
				this.options.collapseBtn(this);
			}

			return this._container;
		},

		onRemove: function() {
		},

		createAlternativesContainer: function() {
			return L.DomUtil.create('div', 'leaflet-routing-alternatives-container');
		},

		setAlternatives: function(routes) {
			var i,
			    alt,
			    altDiv;

			this._clearAlts();

			this._routes = routes;

			for (i = 0; i < this._routes.length; i++) {
				alt = this._routes[i];
				altDiv = this._createAlternative(alt, i);
				this._altContainer.appendChild(altDiv);
				this._altElements.push(altDiv);
			}

			this._selectRoute({route: this._routes[0], alternatives: this._routes.slice(1)});

			return this;
		},

		show: function() {
			L.DomUtil.removeClass(this._container, 'leaflet-routing-container-hide');
		},

		hide: function() {
			L.DomUtil.addClass(this._container, 'leaflet-routing-container-hide');
		},

		_toggle: function() {
			var collapsed = L.DomUtil.hasClass(this._container, 'leaflet-routing-container-hide');
			this[collapsed ? 'show' : 'hide']();
		},

		_createAlternative: function(alt, i) {
			var altDiv = L.DomUtil.create('div', 'leaflet-routing-alt ' +
				this.options.alternativeClassName +
				(i > 0 ? ' leaflet-routing-alt-minimized ' + this.options.minimizedClassName : '')),
				template = this.options.summaryTemplate,
				data = L.extend({
					name: alt.name,
					distance: this._formatter.formatDistance(alt.summary.totalDistance, this.options.totalDistanceRoundingSensitivity),
					time: this._formatter.formatTime(alt.summary.totalTime)
				}, alt);
			altDiv.innerHTML = typeof(template) === 'function' ? template(data) : L.Util.template(template, data);
			L.DomEvent.addListener(altDiv, 'click', this._onAltClicked, this);
			this.on('routeselected', this._selectAlt, this);

			altDiv.appendChild(this._createItineraryContainer(alt));
			return altDiv;
		},

		_clearAlts: function() {
			var el = this._altContainer;
			while (el && el.firstChild) {
				el.removeChild(el.firstChild);
			}

			this._altElements = [];
		},

		_createItineraryContainer: function(r) {
			var container = this._itineraryBuilder.createContainer(),
			    steps = this._itineraryBuilder.createStepsContainer(),
			    i,
			    instr,
			    step,
			    distance,
			    text,
			    icon;

			container.appendChild(steps);

			for (i = 0; i < r.instructions.length; i++) {
				instr = r.instructions[i];
				text = this._formatter.formatInstruction(instr, i);
				distance = this._formatter.formatDistance(instr.distance);
				icon = this._formatter.getIconName(instr, i);
				step = this._itineraryBuilder.createStep(text, distance, icon, steps);

				if(instr.index) {
					this._addRowListeners(step, r.coordinates[instr.index]);
				}
			}

			return container;
		},

		_addRowListeners: function(row, coordinate) {
			L.DomEvent.addListener(row, 'mouseover', function() {
				this._marker = L.circleMarker(coordinate,
					this.options.pointMarkerStyle).addTo(this._map);
			}, this);
			L.DomEvent.addListener(row, 'mouseout', function() {
				if (this._marker) {
					this._map.removeLayer(this._marker);
					delete this._marker;
				}
			}, this);
			L.DomEvent.addListener(row, 'click', function(e) {
				this._map.panTo(coordinate);
				L.DomEvent.stopPropagation(e);
			}, this);
		},

		_onAltClicked: function(e) {
			var altElem = e.target || window.event.srcElement;
			while (!L.DomUtil.hasClass(altElem, 'leaflet-routing-alt')) {
				altElem = altElem.parentElement;
			}

			var j = this._altElements.indexOf(altElem);
			var alts = this._routes.slice();
			var route = alts.splice(j, 1)[0];

			this.fire('routeselected', {
				route: route,
				alternatives: alts
			});
		},

		_selectAlt: function(e) {
			var altElem,
			    j,
			    n,
			    classFn;

			altElem = this._altElements[e.route.routesIndex];

			if (L.DomUtil.hasClass(altElem, 'leaflet-routing-alt-minimized')) {
				for (j = 0; j < this._altElements.length; j++) {
					n = this._altElements[j];
					classFn = j === e.route.routesIndex ? 'removeClass' : 'addClass';
					L.DomUtil[classFn](n, 'leaflet-routing-alt-minimized');
					if (this.options.minimizedClassName) {
						L.DomUtil[classFn](n, this.options.minimizedClassName);
					}

					if (j !== e.route.routesIndex) n.scrollTop = 0;
				}
			}

			L.DomEvent.stop(e);
		},

		_selectRoute: function(routes) {
			if (this._marker) {
				this._map.removeLayer(this._marker);
				delete this._marker;
			}
			this.fire('routeselected', routes);
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./formatter":51,"./itinerary-builder":54}],56:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	module.exports = L.LayerGroup.extend({
		includes: ((typeof L.Evented !== 'undefined' && L.Evented.prototype) || L.Mixin.Events),

		options: {
			styles: [
				{color: 'black', opacity: 0.15, weight: 9}, ///εδω πρεπει να αλλαξω κ'α'τι
				{color: 'white', opacity: 0.8, weight: 6},
				//{color: 'red', opacity: 1, weight: 2}
				{color: 'yellow', opacity: 1, weight: 2}
			],
			missingRouteStyles: [
				{color: 'black', opacity: 0.15, weight: 7},
				{color: 'white', opacity: 0.6, weight: 4},
				{color: 'gray', opacity: 0.8, weight: 2, dashArray: '7,12'}
			],
			addWaypoints: true,
			extendToWaypoints: true,
			missingRouteTolerance: 10
		},

		initialize: function(route, options) {
			L.setOptions(this, options);
			L.LayerGroup.prototype.initialize.call(this, options);
			this._route = route;

			if (this.options.extendToWaypoints) {
				this._extendToWaypoints();
			}

			this._addSegment(
				route.coordinates,
				this.options.styles,
				this.options.addWaypoints);
		},

		getBounds: function() {
			return L.latLngBounds(this._route.coordinates);
		},

		_findWaypointIndices: function() {
			var wps = this._route.inputWaypoints,
			    indices = [],
			    i;
			for (i = 0; i < wps.length; i++) {
				indices.push(this._findClosestRoutePoint(wps[i].latLng));
			}

			return indices;
		},

		_findClosestRoutePoint: function(latlng) {
			var minDist = Number.MAX_VALUE,
				minIndex,
			    i,
			    d;

			for (i = this._route.coordinates.length - 1; i >= 0 ; i--) {
				// TODO: maybe do this in pixel space instead?
				d = latlng.distanceTo(this._route.coordinates[i]);
				if (d < minDist) {
					minIndex = i;
					minDist = d;
				}
			}

			return minIndex;
		},

		_extendToWaypoints: function() {
			var wps = this._route.inputWaypoints,
				wpIndices = this._getWaypointIndices(),
			    i,
			    wpLatLng,
			    routeCoord;

			for (i = 0; i < wps.length; i++) {
				wpLatLng = wps[i].latLng;
				routeCoord = L.latLng(this._route.coordinates[wpIndices[i]]);
				if (wpLatLng.distanceTo(routeCoord) >
					this.options.missingRouteTolerance) {
					this._addSegment([wpLatLng, routeCoord],
						this.options.missingRouteStyles);
				}
			}
		},

		_addSegment: function(coords, styles, mouselistener) {
			var i,
				pl;

			for (i = 0; i < styles.length; i++) {
				pl = L.polyline(coords, styles[i]);
				this.addLayer(pl);
				if (mouselistener) {
					pl.on('mousedown', this._onLineTouched, this);
				}
			}
		},

		_findNearestWpBefore: function(i) {
			var wpIndices = this._getWaypointIndices(),
				j = wpIndices.length - 1;
			while (j >= 0 && wpIndices[j] > i) {
				j--;
			}

			return j;
		},

		_onLineTouched: function(e) {
			var afterIndex = this._findNearestWpBefore(this._findClosestRoutePoint(e.latlng));
			this.fire('linetouched', {
				afterIndex: afterIndex,
				latlng: e.latlng
			});
			L.DomEvent.stop(e);
		},

		_getWaypointIndices: function() {
			if (!this._wpIndices) {
				this._wpIndices = this._route.waypointIndices || this._findWaypointIndices();
			}

			return this._wpIndices;
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}],57:[function(_dereq_,module,exports){
/* 
   NOTICE
   Since version 3.2.5, the functionality in this file is by
   default NOT used for localizing OSRM instructions.
   Instead, we rely on the module osrm-text-instructions (https://github.com/Project-OSRM/osrm-text-instructions/).
   
   This file can still be used for other routing backends, or if you specify the
   stepToText option in the OSRMv1 class.
*/

(function() {
	//'use strict';

	var spanish = {
		directions: {
			N: 'norte',
			NE: 'noreste',
			E: 'este',
			SE: 'sureste',
			S: 'sur',
			SW: 'suroeste',
			W: 'oeste',
			NW: 'noroeste',
			SlightRight: 'leve giro a la derecha',
			Right: 'derecha',
			SharpRight: 'giro pronunciado a la derecha',
			SlightLeft: 'leve giro a la izquierda',
			Left: 'izquierda',
			SharpLeft: 'giro pronunciado a la izquierda',
			Uturn: 'media vuelta'
		},
		instructions: {
			// instruction, postfix if the road is named
			'Head':
				['Derecho {dir}', ' sobre {road}'],
			'Continue':
				['Continuar {dir}', ' en {road}'],
			'TurnAround':
				['Dar vuelta'],
			'WaypointReached':
				['Llegó a un punto del camino'],
			'Roundabout':
				['Tomar {exitStr} salida en la rotonda', ' en {road}'],
			'DestinationReached':
				['Llegada a destino'],
			'Fork': ['En el cruce gira a {modifier}', ' hacia {road}'],
			'Merge': ['Incorpórate {modifier}', ' hacia {road}'],
			'OnRamp': ['Gira {modifier} en la salida', ' hacia {road}'],
			'OffRamp': ['Toma la salida {modifier}', ' hacia {road}'],
			'EndOfRoad': ['Gira {modifier} al final de la carretera', ' hacia {road}'],
			'Onto': 'hacia {road}'
		},
		formatOrder: function(n) {
			return n + 'º';
		},
		ui: {
			startPlaceholder: 'Inicio',
			viaPlaceholder: 'Via {viaNumber}',
			endPlaceholder: 'Destino'
		},
		units: {
			meters: 'm',
			kilometers: 'km',
			yards: 'yd',
			miles: 'mi',
			hours: 'h',
			minutes: 'min',
			seconds: 's'
		}
	};

	L.Routing = L.Routing || {};

	var Localization = L.Class.extend({
		initialize: function(langs) {
			this._langs = L.Util.isArray(langs) ? langs.slice() : [langs, 'en'];

			for (var i = 0, l = this._langs.length; i < l; i++) {
				var generalizedCode = /([A-Za-z]+)/.exec(this._langs[i])[1]
				if (!Localization[this._langs[i]]) {
					if (Localization[generalizedCode]) {
						this._langs[i] = generalizedCode;
					} else {
						throw new Error('No localization for language "' + this._langs[i] + '".');
					}
				}
			}
		},

		localize: function(keys) {
			var dict,
				key,
				value;

			keys = L.Util.isArray(keys) ? keys : [keys];

			for (var i = 0, l = this._langs.length; i < l; i++) {
				dict = Localization[this._langs[i]];
				for (var j = 0, nKeys = keys.length; dict && j < nKeys; j++) {
					key = keys[j];
					value = dict[key];
					dict = value;
				}

				if (value) {
					return value;
				}
			}
		}
	});

	module.exports = L.extend(Localization, {
		'en': {
			directions: {
				N: 'north',
				NE: 'northeast',
				E: 'east',
				SE: 'southeast',
				S: 'south',
				SW: 'southwest',
				W: 'west',
				NW: 'northwest',
				SlightRight: 'slight right',
				Right: 'right',
				SharpRight: 'sharp right',
				SlightLeft: 'slight left',
				Left: 'left',
				SharpLeft: 'sharp left',
				Uturn: 'Turn around'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Head {dir}', ' on {road}'],
				'Continue':
					['Continue {dir}'],
				'TurnAround':
					['Turn around'],
				'WaypointReached':
					['Waypoint reached'],
				'Roundabout':
					['Take the {exitStr} exit in the roundabout', ' onto {road}'],
				'DestinationReached':
					['Destination reached'],
				'Fork': ['At the fork, turn {modifier}', ' onto {road}'],
				'Merge': ['Merge {modifier}', ' onto {road}'],
				'OnRamp': ['Turn {modifier} on the ramp', ' onto {road}'],
				'OffRamp': ['Take the ramp on the {modifier}', ' onto {road}'],
				'EndOfRoad': ['Turn {modifier} at the end of the road', ' onto {road}'],
				'Onto': 'onto {road}'
			},
			formatOrder: function(n) {
				var i = n % 10 - 1,
				suffix = ['st', 'nd', 'rd'];

				return suffix[i] ? n + suffix[i] : n + 'th';
			},
			ui: {
				startPlaceholder: 'Start',
				viaPlaceholder: 'Via {viaNumber}',
				endPlaceholder: 'End'
			},
			units: {
				meters: 'm',
				kilometers: 'km',
				yards: 'yd',
				miles: 'mi',
				hours: 'h',
				minutes: 'min',
				seconds: 's'
			}
		},

		'de': {
			directions: {
				N: 'Norden',
				NE: 'Nordosten',
				E: 'Osten',
				SE: 'Südosten',
				S: 'Süden',
				SW: 'Südwesten',
				W: 'Westen',
				NW: 'Nordwesten',
				SlightRight: 'leicht rechts',
				Right: 'rechts',
				SharpRight: 'scharf rechts',
				SlightLeft: 'leicht links',
				Left: 'links',
				SharpLeft: 'scharf links',
				Uturn: 'Wenden'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Richtung {dir}', ' auf {road}'],
				'Continue':
					['Geradeaus Richtung {dir}', ' auf {road}'],
				'SlightRight':
					['Leicht rechts abbiegen', ' auf {road}'],
				'Right':
					['Rechts abbiegen', ' auf {road}'],
				'SharpRight':
					['Scharf rechts abbiegen', ' auf {road}'],
				'TurnAround':
					['Wenden'],
				'SharpLeft':
					['Scharf links abbiegen', ' auf {road}'],
				'Left':
					['Links abbiegen', ' auf {road}'],
				'SlightLeft':
					['Leicht links abbiegen', ' auf {road}'],
				'WaypointReached':
					['Zwischenhalt erreicht'],
				'Roundabout':
					['Nehmen Sie die {exitStr} Ausfahrt im Kreisverkehr', ' auf {road}'],
				'DestinationReached':
					['Sie haben ihr Ziel erreicht'],
				'Fork': ['An der Kreuzung {modifier}', ' auf {road}'],
				'Merge': ['Fahren Sie {modifier} weiter', ' auf {road}'],
				'OnRamp': ['Fahren Sie {modifier} auf die Auffahrt', ' auf {road}'],
				'OffRamp': ['Nehmen Sie die Ausfahrt {modifier}', ' auf {road}'],
				'EndOfRoad': ['Fahren Sie {modifier} am Ende der Straße', ' auf {road}'],
				'Onto': 'auf {road}'
			},
			formatOrder: function(n) {
				return n + '.';
			},
			ui: {
				startPlaceholder: 'Start',
				viaPlaceholder: 'Via {viaNumber}',
				endPlaceholder: 'Ziel'
			}
		},

		'sv': {
			directions: {
				N: 'norr',
				NE: 'nordost',
				E: 'öst',
				SE: 'sydost',
				S: 'syd',
				SW: 'sydväst',
				W: 'väst',
				NW: 'nordväst',
				SlightRight: 'svagt höger',
				Right: 'höger',
				SharpRight: 'skarpt höger',
				SlightLeft: 'svagt vänster',
				Left: 'vänster',
				SharpLeft: 'skarpt vänster',
				Uturn: 'Vänd'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Åk åt {dir}', ' till {road}'],
				'Continue':
					['Fortsätt {dir}'],
				'SlightRight':
					['Svagt höger', ' till {road}'],
				'Right':
					['Sväng höger', ' till {road}'],
				'SharpRight':
					['Skarpt höger', ' till {road}'],
				'TurnAround':
					['Vänd'],
				'SharpLeft':
					['Skarpt vänster', ' till {road}'],
				'Left':
					['Sväng vänster', ' till {road}'],
				'SlightLeft':
					['Svagt vänster', ' till {road}'],
				'WaypointReached':
					['Viapunkt nådd'],
				'Roundabout':
					['Tag {exitStr} avfarten i rondellen', ' till {road}'],
				'DestinationReached':
					['Framme vid resans mål'],
				'Fork': ['Tag av {modifier}', ' till {road}'],
				'Merge': ['Anslut {modifier} ', ' till {road}'],
				'OnRamp': ['Tag påfarten {modifier}', ' till {road}'],
				'OffRamp': ['Tag avfarten {modifier}', ' till {road}'],
				'EndOfRoad': ['Sväng {modifier} vid vägens slut', ' till {road}'],
				'Onto': 'till {road}'
			},
			formatOrder: function(n) {
				return ['första', 'andra', 'tredje', 'fjärde', 'femte',
					'sjätte', 'sjunde', 'åttonde', 'nionde', 'tionde'
					/* Can't possibly be more than ten exits, can there? */][n - 1];
			},
			ui: {
				startPlaceholder: 'Från',
				viaPlaceholder: 'Via {viaNumber}',
				endPlaceholder: 'Till'
			}
		},

		'es': spanish,
		'sp': spanish,
		
		'nl': {
			directions: {
				N: 'noordelijke',
				NE: 'noordoostelijke',
				E: 'oostelijke',
				SE: 'zuidoostelijke',
				S: 'zuidelijke',
				SW: 'zuidewestelijke',
				W: 'westelijke',
				NW: 'noordwestelijke'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Vertrek in {dir} richting', ' de {road} op'],
				'Continue':
					['Ga in {dir} richting', ' de {road} op'],
				'SlightRight':
					['Volg de weg naar rechts', ' de {road} op'],
				'Right':
					['Ga rechtsaf', ' de {road} op'],
				'SharpRight':
					['Ga scherpe bocht naar rechts', ' de {road} op'],
				'TurnAround':
					['Keer om'],
				'SharpLeft':
					['Ga scherpe bocht naar links', ' de {road} op'],
				'Left':
					['Ga linksaf', ' de {road} op'],
				'SlightLeft':
					['Volg de weg naar links', ' de {road} op'],
				'WaypointReached':
					['Aangekomen bij tussenpunt'],
				'Roundabout':
					['Neem de {exitStr} afslag op de rotonde', ' de {road} op'],
				'DestinationReached':
					['Aangekomen op eindpunt'],
			},
			formatOrder: function(n) {
				if (n === 1 || n >= 20) {
					return n + 'ste';
				} else {
					return n + 'de';
				}
			},
			ui: {
				startPlaceholder: 'Vertrekpunt',
				viaPlaceholder: 'Via {viaNumber}',
				endPlaceholder: 'Bestemming'
			}
		},
		'fr': {
			directions: {
				N: 'nord',
				NE: 'nord-est',
				E: 'est',
				SE: 'sud-est',
				S: 'sud',
				SW: 'sud-ouest',
				W: 'ouest',
				NW: 'nord-ouest'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Tout droit au {dir}', ' sur {road}'],
				'Continue':
					['Continuer au {dir}', ' sur {road}'],
				'SlightRight':
					['Légèrement à droite', ' sur {road}'],
				'Right':
					['A droite', ' sur {road}'],
				'SharpRight':
					['Complètement à droite', ' sur {road}'],
				'TurnAround':
					['Faire demi-tour'],
				'SharpLeft':
					['Complètement à gauche', ' sur {road}'],
				'Left':
					['A gauche', ' sur {road}'],
				'SlightLeft':
					['Légèrement à gauche', ' sur {road}'],
				'WaypointReached':
					['Point d\'étape atteint'],
				'Roundabout':
					['Au rond-point, prenez la {exitStr} sortie', ' sur {road}'],
				'DestinationReached':
					['Destination atteinte'],
			},
			formatOrder: function(n) {
				return n + 'º';
			},
			ui: {
				startPlaceholder: 'Départ',
				viaPlaceholder: 'Intermédiaire {viaNumber}',
				endPlaceholder: 'Arrivée'
			}
		},
		'it': {
			directions: {
				N: 'nord',
				NE: 'nord-est',
				E: 'est',
				SE: 'sud-est',
				S: 'sud',
				SW: 'sud-ovest',
				W: 'ovest',
				NW: 'nord-ovest'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Dritto verso {dir}', ' su {road}'],
				'Continue':
					['Continuare verso {dir}', ' su {road}'],
				'SlightRight':
					['Mantenere la destra', ' su {road}'],
				'Right':
					['A destra', ' su {road}'],
				'SharpRight':
					['Strettamente a destra', ' su {road}'],
				'TurnAround':
					['Fare inversione di marcia'],
				'SharpLeft':
					['Strettamente a sinistra', ' su {road}'],
				'Left':
					['A sinistra', ' sur {road}'],
				'SlightLeft':
					['Mantenere la sinistra', ' su {road}'],
				'WaypointReached':
					['Punto di passaggio raggiunto'],
				'Roundabout':
					['Alla rotonda, prendere la {exitStr} uscita'],
				'DestinationReached':
					['Destinazione raggiunta'],
			},
			formatOrder: function(n) {
				return n + 'º';
			},
			ui: {
				startPlaceholder: 'Partenza',
				viaPlaceholder: 'Intermedia {viaNumber}',
				endPlaceholder: 'Destinazione'
			}
		},
		'pt': {
			directions: {
				N: 'norte',
				NE: 'nordeste',
				E: 'leste',
				SE: 'sudeste',
				S: 'sul',
				SW: 'sudoeste',
				W: 'oeste',
				NW: 'noroeste',
				SlightRight: 'curva ligeira a direita',
				Right: 'direita',
				SharpRight: 'curva fechada a direita',
				SlightLeft: 'ligeira a esquerda',
				Left: 'esquerda',
				SharpLeft: 'curva fechada a esquerda',
				Uturn: 'Meia volta'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Siga {dir}', ' na {road}'],
				'Continue':
					['Continue {dir}', ' na {road}'],
				'SlightRight':
					['Curva ligeira a direita', ' na {road}'],
				'Right':
					['Curva a direita', ' na {road}'],
				'SharpRight':
					['Curva fechada a direita', ' na {road}'],
				'TurnAround':
					['Retorne'],
				'SharpLeft':
					['Curva fechada a esquerda', ' na {road}'],
				'Left':
					['Curva a esquerda', ' na {road}'],
				'SlightLeft':
					['Curva ligueira a esquerda', ' na {road}'],
				'WaypointReached':
					['Ponto de interesse atingido'],
				'Roundabout':
					['Pegue a {exitStr} saída na rotatória', ' na {road}'],
				'DestinationReached':
					['Destino atingido'],
				'Fork': ['Na encruzilhada, vire a {modifier}', ' na {road}'],
				'Merge': ['Entre à {modifier}', ' na {road}'],
				'OnRamp': ['Vire {modifier} na rampa', ' na {road}'],
				'OffRamp': ['Entre na rampa na {modifier}', ' na {road}'],
				'EndOfRoad': ['Vire {modifier} no fim da rua', ' na {road}'],
				'Onto': 'na {road}'
			},
			formatOrder: function(n) {
				return n + 'º';
			},
			ui: {
				startPlaceholder: 'Origem',
				viaPlaceholder: 'Intermédio {viaNumber}',
				endPlaceholder: 'Destino'
			}
		},
		'sk': {
			directions: {
				N: 'sever',
				NE: 'serverovýchod',
				E: 'východ',
				SE: 'juhovýchod',
				S: 'juh',
				SW: 'juhozápad',
				W: 'západ',
				NW: 'serverozápad'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Mierte na {dir}', ' na {road}'],
				'Continue':
					['Pokračujte na {dir}', ' na {road}'],
				'SlightRight':
					['Mierne doprava', ' na {road}'],
				'Right':
					['Doprava', ' na {road}'],
				'SharpRight':
					['Prudko doprava', ' na {road}'],
				'TurnAround':
					['Otočte sa'],
				'SharpLeft':
					['Prudko doľava', ' na {road}'],
				'Left':
					['Doľava', ' na {road}'],
				'SlightLeft':
					['Mierne doľava', ' na {road}'],
				'WaypointReached':
					['Ste v prejazdovom bode.'],
				'Roundabout':
					['Odbočte na {exitStr} výjazde', ' na {road}'],
				'DestinationReached':
					['Prišli ste do cieľa.'],
			},
			formatOrder: function(n) {
				var i = n % 10 - 1,
				suffix = ['.', '.', '.'];

				return suffix[i] ? n + suffix[i] : n + '.';
			},
			ui: {
				startPlaceholder: 'Začiatok',
				viaPlaceholder: 'Cez {viaNumber}',
				endPlaceholder: 'Koniec'
			}
		},
		'el': {
			directions: {
				N: 'βόρεια',
				NE: 'βορειοανατολικά',
				E: 'ανατολικά',
				SE: 'νοτιοανατολικά',
				S: 'νότια',
				SW: 'νοτιοδυτικά',
				W: 'δυτικά',
				NW: 'βορειοδυτικά'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Κατευθυνθείτε {dir}', ' στην {road}'],
				'Continue':
					['Συνεχίστε {dir}', ' στην {road}'],
				'SlightRight':
					['Ελαφρώς δεξιά', ' στην {road}'],
				'Right':
					['Δεξιά', ' στην {road}'],
				'SharpRight':
					['Απότομη δεξιά στροφή', ' στην {road}'],
				'TurnAround':
					['Κάντε αναστροφή'],
				'SharpLeft':
					['Απότομη αριστερή στροφή', ' στην {road}'],
				'Left':
					['Αριστερά', ' στην {road}'],
				'SlightLeft':
					['Ελαφρώς αριστερά', ' στην {road}'],
				'WaypointReached':
					['Φτάσατε στο σημείο αναφοράς'],
				'Roundabout':
					['Ακολουθήστε την {exitStr} έξοδο στο κυκλικό κόμβο', ' στην {road}'],
				'DestinationReached':
					['Φτάσατε στον προορισμό σας'],
			},
			formatOrder: function(n) {
				return n + 'º';
			},
			ui: {
				startPlaceholder: 'Αφετηρία',
				viaPlaceholder: 'μέσω {viaNumber}',
				endPlaceholder: 'Προορισμός'
			}
		},
		'ca': {
			directions: {
				N: 'nord',
				NE: 'nord-est',
				E: 'est',
				SE: 'sud-est',
				S: 'sud',
				SW: 'sud-oest',
				W: 'oest',
				NW: 'nord-oest',
				SlightRight: 'lleu gir a la dreta',
				Right: 'dreta',
				SharpRight: 'gir pronunciat a la dreta',
				SlightLeft: 'gir pronunciat a l\'esquerra',
				Left: 'esquerra',
				SharpLeft: 'lleu gir a l\'esquerra',
				Uturn: 'mitja volta'
			},
			instructions: {
				'Head':
					['Recte {dir}', ' sobre {road}'],
				'Continue':
					['Continuar {dir}'],
				'TurnAround':
					['Donar la volta'],
				'WaypointReached':
					['Ha arribat a un punt del camí'],
				'Roundabout':
					['Agafar {exitStr} sortida a la rotonda', ' a {road}'],
				'DestinationReached':
					['Arribada al destí'],
				'Fork': ['A la cruïlla gira a la {modifier}', ' cap a {road}'],
				'Merge': ['Incorpora\'t {modifier}', ' a {road}'],
				'OnRamp': ['Gira {modifier} a la sortida', ' cap a {road}'],
				'OffRamp': ['Pren la sortida {modifier}', ' cap a {road}'],
				'EndOfRoad': ['Gira {modifier} al final de la carretera', ' cap a {road}'],
				'Onto': 'cap a {road}'
			},
			formatOrder: function(n) {
				return n + 'º';
			},
			ui: {
				startPlaceholder: 'Origen',
				viaPlaceholder: 'Via {viaNumber}',
				endPlaceholder: 'Destí'
			},
			units: {
				meters: 'm',
				kilometers: 'km',
				yards: 'yd',
				miles: 'mi',
				hours: 'h',
				minutes: 'min',
				seconds: 's'
			}
		},
		'ru': {
			directions: {
				N: 'север',
				NE: 'северовосток',
				E: 'восток',
				SE: 'юговосток',
				S: 'юг',
				SW: 'югозапад',
				W: 'запад',
				NW: 'северозапад',
				SlightRight: 'плавно направо',
				Right: 'направо',
				SharpRight: 'резко направо',
				SlightLeft: 'плавно налево',
				Left: 'налево',
				SharpLeft: 'резко налево',
				Uturn: 'развернуться'
			},
			instructions: {
				'Head':
					['Начать движение на {dir}', ' по {road}'],
				'Continue':
					['Продолжать движение на {dir}', ' по {road}'],
				'SlightRight':
					['Плавный поворот направо', ' на {road}'],
				'Right':
					['Направо', ' на {road}'],
				'SharpRight':
					['Резкий поворот направо', ' на {road}'],
				'TurnAround':
					['Развернуться'],
				'SharpLeft':
					['Резкий поворот налево', ' на {road}'],
				'Left':
					['Поворот налево', ' на {road}'],
				'SlightLeft':
					['Плавный поворот налево', ' на {road}'],
				'WaypointReached':
					['Точка достигнута'],
				'Roundabout':
					['{exitStr} съезд с кольца', ' на {road}'],
				'DestinationReached':
					['Окончание маршрута'],
				'Fork': ['На развилке поверните {modifier}', ' на {road}'],
				'Merge': ['Перестройтесь {modifier}', ' на {road}'],
				'OnRamp': ['Поверните {modifier} на съезд', ' на {road}'],
				'OffRamp': ['Съезжайте на {modifier}', ' на {road}'],
				'EndOfRoad': ['Поверните {modifier} в конце дороги', ' на {road}'],
				'Onto': 'на {road}'
			},
			formatOrder: function(n) {
				return n + '-й';
			},
			ui: {
				startPlaceholder: 'Начало',
				viaPlaceholder: 'Через {viaNumber}',
				endPlaceholder: 'Конец'
			},
			units: {
				meters: 'м',
				kilometers: 'км',
				yards: 'ярд',
				miles: 'ми',
				hours: 'ч',
				minutes: 'м',
				seconds: 'с'
			}
		},
                
                'pl': {
			directions: {
				N: 'północ',
				NE: 'północny wschód',
				E: 'wschód',
				SE: 'południowy wschód',
				S: 'południe',
				SW: 'południowy zachód',
				W: 'zachód',
				NW: 'północny zachód',
				SlightRight: 'lekko w prawo',
				Right: 'w prawo',
				SharpRight: 'ostro w prawo',
				SlightLeft: 'lekko w lewo',
				Left: 'w lewo',
				SharpLeft: 'ostro w lewo',
				Uturn: 'zawróć'
			},
			instructions: {
				// instruction, postfix if the road is named
				'Head':
					['Kieruj się na {dir}', ' na {road}'],
				'Continue':
					['Jedź dalej przez {dir}'],
				'TurnAround':
					['Zawróć'],
				'WaypointReached':
					['Punkt pośredni'],
				'Roundabout':
					['Wyjedź {exitStr} zjazdem na rondzie', ' na {road}'],
				'DestinationReached':
					['Dojechano do miejsca docelowego'],
				'Fork': ['Na rozwidleniu {modifier}', ' na {road}'],
				'Merge': ['Zjedź {modifier}', ' na {road}'],
				'OnRamp': ['Wjazd {modifier}', ' na {road}'],
				'OffRamp': ['Zjazd {modifier}', ' na {road}'],
				'EndOfRoad': ['Skręć {modifier} na końcu drogi', ' na {road}'],
				'Onto': 'na {road}'
			},
			formatOrder: function(n) {
				return n + '.';
			},
			ui: {
				startPlaceholder: 'Początek',
				viaPlaceholder: 'Przez {viaNumber}',
				endPlaceholder: 'Koniec'
			},
			units: {
				meters: 'm',
				kilometers: 'km',
				yards: 'yd',
				miles: 'mi',
				hours: 'godz',
				minutes: 'min',
				seconds: 's'
			}
		}
	});
})();

},{}],58:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	var OSRMv1 = _dereq_('./osrm-v1');

	/**
	 * Works against OSRM's new API in version 5.0; this has
	 * the API version v1.
	 */
	module.exports = OSRMv1.extend({
		options: {
			serviceUrl: 'https://api.mapbox.com/directions/v5',
			profile: 'mapbox/driving',
			useHints: false
		},

		initialize: function(accessToken, options) {
			L.Routing.OSRMv1.prototype.initialize.call(this, options);
			this.options.requestParameters = this.options.requestParameters || {};
			/* jshint camelcase: false */
			this.options.requestParameters.access_token = accessToken;
			/* jshint camelcase: true */
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./osrm-v1":59}],59:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null),
		corslite = _dereq_('@mapbox/corslite'),
		polyline = _dereq_('@mapbox/polyline'),
		osrmTextInstructions = _dereq_('osrm-text-instructions')('v5');

	// Ignore camelcase naming for this file, since OSRM's API uses
	// underscores.
	/* jshint camelcase: false */

	var Waypoint = _dereq_('./waypoint');

	/**
	 * Works against OSRM's new API in version 5.0; this has
	 * the API version v1.
	 */
	module.exports = L.Class.extend({
		options: {
			serviceUrl: 'https://router.project-osrm.org/route/v1',
			profile: 'driving',
			timeout: 30 * 1000,
			routingOptions: {
				alternatives: true,
				steps: true
			},
			polylinePrecision: 5,
			useHints: true,
			suppressDemoServerWarning: false,
			language: 'en'
		},

		initialize: function(options) {
			L.Util.setOptions(this, options);
			this._hints = {
				locations: {}
			};

			if (!this.options.suppressDemoServerWarning &&
				this.options.serviceUrl.indexOf('//router.project-osrm.org') >= 0) {
				console.warn('You are using OSRM\'s demo server. ' +
					'Please note that it is **NOT SUITABLE FOR PRODUCTION USE**.\n' +
					'Refer to the demo server\'s usage policy: ' +
					'https://github.com/Project-OSRM/osrm-backend/wiki/Api-usage-policy\n\n' +
					'To change, set the serviceUrl option.\n\n' +
					'Please do not report issues with this server to neither ' +
					'Leaflet Routing Machine or OSRM - it\'s for\n' +
					'demo only, and will sometimes not be available, or work in ' +
					'unexpected ways.\n\n' +
					'Please set up your own OSRM server, or use a paid service ' +
					'provider for production.');
			}
		},

		route: function(waypoints, callback, context, options) {
			var timedOut = false,
				wps = [],
				url,
				timer,
				wp,
				i,
				xhr;

			options = L.extend({}, this.options.routingOptions, options);
			url = this.buildRouteUrl(waypoints, options);
			if (this.options.requestParameters) {
				url += L.Util.getParamString(this.options.requestParameters, url);
			}

			timer = setTimeout(function() {
				timedOut = true;
				callback.call(context || callback, {
					status: -1,
					message: 'OSRM request timed out.'
				});
			}, this.options.timeout);

			// Create a copy of the waypoints, since they
			// might otherwise be asynchronously modified while
			// the request is being processed.
			for (i = 0; i < waypoints.length; i++) {
				wp = waypoints[i];
				wps.push(new Waypoint(wp.latLng, wp.name, wp.options));
			}

			return xhr = corslite(url, L.bind(function(err, resp) {
				var data,
					error =  {};

				clearTimeout(timer);
				if (!timedOut) {
					if (!err) {
						try {
							data = JSON.parse(resp.responseText);
							try {
								return this._routeDone(data, wps, options, callback, context);
							} catch (ex) {
								error.status = -3;
								error.message = ex.toString();
							}
						} catch (ex) {
							error.status = -2;
							error.message = 'Error parsing OSRM response: ' + ex.toString();
						}
					} else {
						error.message = 'HTTP request failed: ' + err.type +
							(err.target && err.target.status ? ' HTTP ' + err.target.status + ': ' + err.target.statusText : '');
						error.url = url;
						error.status = -1;
						error.target = err;
					}

					callback.call(context || callback, error);
				} else {
					xhr.abort();
				}
			}, this));
		},

		requiresMoreDetail: function(route, zoom, bounds) {
			if (!route.properties.isSimplified) {
				return false;
			}

			var waypoints = route.inputWaypoints,
				i;
			for (i = 0; i < waypoints.length; ++i) {
				if (!bounds.contains(waypoints[i].latLng)) {
					return true;
				}
			}

			return false;
		},

		_routeDone: function(response, inputWaypoints, options, callback, context) {
			var alts = [],
			    actualWaypoints,
			    i,
			    route;

			context = context || callback;
			if (response.code !== 'Ok') {
				callback.call(context, {
					status: response.code
				});
				return;
			}

			actualWaypoints = this._toWaypoints(inputWaypoints, response.waypoints);

			for (i = 0; i < response.routes.length; i++) {
				route = this._convertRoute(response.routes[i]);
				route.inputWaypoints = inputWaypoints;
				route.waypoints = actualWaypoints;
				route.properties = {isSimplified: !options || !options.geometryOnly || options.simplifyGeometry};
				alts.push(route);
			}

			this._saveHintData(response.waypoints, inputWaypoints);

			callback.call(context, null, alts);
		},

		_convertRoute: function(responseRoute) {
			var result = {
					name: '',
					coordinates: [],
					instructions: [],
					summary: {
						totalDistance: responseRoute.distance,
						totalTime: responseRoute.duration
					}
				},
				legNames = [],
				waypointIndices = [],
				index = 0,
				legCount = responseRoute.legs.length,
				hasSteps = responseRoute.legs[0].steps.length > 0,
				i,
				j,
				leg,
				step,
				geometry,
				type,
				modifier,
				text,
				stepToText;

			if (this.options.stepToText) {
				stepToText = this.options.stepToText;
			} else {
				stepToText = L.bind(osrmTextInstructions.compile, osrmTextInstructions, this.options.language);
			}

			for (i = 0; i < legCount; i++) {
				leg = responseRoute.legs[i];
				legNames.push(leg.summary && leg.summary.charAt(0).toUpperCase() + leg.summary.substring(1));
				for (j = 0; j < leg.steps.length; j++) {
					step = leg.steps[j];
					geometry = this._decodePolyline(step.geometry);
					result.coordinates.push.apply(result.coordinates, geometry);
					type = this._maneuverToInstructionType(step.maneuver, i === legCount - 1);
					modifier = this._maneuverToModifier(step.maneuver);
					text = stepToText(step, {legCount: legCount, legIndex: i});

					if (type) {
						if ((i == 0 && step.maneuver.type == 'depart') || step.maneuver.type == 'arrive') {
							waypointIndices.push(index);
						}

						result.instructions.push({
							type: type,
							distance: step.distance,
							time: step.duration,
							road: step.name,
							direction: this._bearingToDirection(step.maneuver.bearing_after),
							exit: step.maneuver.exit,
							index: index,
							mode: step.mode,
							modifier: modifier,
							text: text
						});
					}

					index += geometry.length;
				}
			}

			result.name = legNames.join(', ');
			if (!hasSteps) {
				result.coordinates = this._decodePolyline(responseRoute.geometry);
			} else {
				result.waypointIndices = waypointIndices;
			}

			return result;
		},

		_bearingToDirection: function(bearing) {
			var oct = Math.round(bearing / 45) % 8;
			return ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'][oct];
		},

		_maneuverToInstructionType: function(maneuver, lastLeg) {
			switch (maneuver.type) {
			case 'new name':
				return 'Continue';
			case 'depart':
				return 'Head';
			case 'arrive':
				return lastLeg ? 'DestinationReached' : 'WaypointReached';
			case 'roundabout':
			case 'rotary':
				return 'Roundabout';
			case 'merge':
			case 'fork':
			case 'on ramp':
			case 'off ramp':
			case 'end of road':
				return this._camelCase(maneuver.type);
			// These are all reduced to the same instruction in the current model
			//case 'turn':
			//case 'ramp': // deprecated in v5.1
			default:
				return this._camelCase(maneuver.modifier);
			}
		},

		_maneuverToModifier: function(maneuver) {
			var modifier = maneuver.modifier;

			switch (maneuver.type) {
			case 'merge':
			case 'fork':
			case 'on ramp':
			case 'off ramp':
			case 'end of road':
				modifier = this._leftOrRight(modifier);
			}

			return modifier && this._camelCase(modifier);
		},

		_camelCase: function(s) {
			var words = s.split(' '),
				result = '';
			for (var i = 0, l = words.length; i < l; i++) {
				result += words[i].charAt(0).toUpperCase() + words[i].substring(1);
			}

			return result;
		},

		_leftOrRight: function(d) {
			return d.indexOf('left') >= 0 ? 'Left' : 'Right';
		},

		_decodePolyline: function(routeGeometry) {
			var cs = polyline.decode(routeGeometry, this.options.polylinePrecision),
				result = new Array(cs.length),
				i;
			for (i = cs.length - 1; i >= 0; i--) {
				result[i] = L.latLng(cs[i]);
			}

			return result;
		},

		_toWaypoints: function(inputWaypoints, vias) {
			var wps = [],
			    i,
			    viaLoc;
			for (i = 0; i < vias.length; i++) {
				viaLoc = vias[i].location;
				wps.push(new Waypoint(L.latLng(viaLoc[1], viaLoc[0]),
				                            inputWaypoints[i].name,
											inputWaypoints[i].options));
			}

			return wps;
		},

		buildRouteUrl: function(waypoints, options) {
			var locs = [],
				hints = [],
				wp,
				latLng,
			    computeInstructions,
			    computeAlternative = true;

			for (var i = 0; i < waypoints.length; i++) {
				wp = waypoints[i];
				latLng = wp.latLng;
				locs.push(latLng.lng + ',' + latLng.lat);
				hints.push(this._hints.locations[this._locationKey(latLng)] || '');
			}

			computeInstructions =
				true;

			return this.options.serviceUrl + '/' + this.options.profile + '/' +
				locs.join(';') + '?' +
				(options.geometryOnly ? (options.simplifyGeometry ? '' : 'overview=full') : 'overview=false') +
				'&alternatives=' + computeAlternative.toString() +
				'&steps=' + computeInstructions.toString() +
				(this.options.useHints ? '&hints=' + hints.join(';') : '') +
				(options.allowUTurns ? '&continue_straight=' + !options.allowUTurns : '');
		},

		_locationKey: function(location) {
			return location.lat + ',' + location.lng;
		},

		_saveHintData: function(actualWaypoints, waypoints) {
			var loc;
			this._hints = {
				locations: {}
			};
			for (var i = actualWaypoints.length - 1; i >= 0; i--) {
				loc = waypoints[i].latLng;
				this._hints.locations[this._locationKey(loc)] = actualWaypoints[i].hint;
			}
		},
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./waypoint":61,"@mapbox/corslite":1,"@mapbox/polyline":2,"osrm-text-instructions":3}],60:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);
	var GeocoderElement = _dereq_('./geocoder-element');
	var Waypoint = _dereq_('./waypoint');

	module.exports = (L.Layer || L.Class).extend({
		includes: ((typeof L.Evented !== 'undefined' && L.Evented.prototype) || L.Mixin.Events),

		options: {
			dragStyles: [
				{color: 'black', opacity: 0.15, weight: 9},
				{color: 'white', opacity: 0.8, weight: 6},
				{color: 'red', opacity: 1, weight: 2, dashArray: '7,12'}
			],
			draggableWaypoints: true,
			routeWhileDragging: false,
			addWaypoints: true,
			reverseWaypoints: false,
			addButtonClassName: '',
			language: 'en',
			createGeocoderElement: function(wp, i, nWps, plan) {
				return new GeocoderElement(wp, i, nWps, plan);
			},
			createMarker: function(i, wp) {
				var options = {
						draggable: this.draggableWaypoints
					},
				    marker = L.marker(wp.latLng, options);

				return marker;
			},
			geocodersClassName: ''
		},

		initialize: function(waypoints, options) {
			L.Util.setOptions(this, options);
			this._waypoints = [];
			this.setWaypoints(waypoints);
		},

		isReady: function() {
			var i;
			for (i = 0; i < this._waypoints.length; i++) {
				if (!this._waypoints[i].latLng) {
					return false;
				}
			}

			return true;
		},

		getWaypoints: function() {
			var i,
				wps = [];

			for (i = 0; i < this._waypoints.length; i++) {
				wps.push(this._waypoints[i]);
			}

			return wps;
		},

		setWaypoints: function(waypoints) {
			var args = [0, this._waypoints.length].concat(waypoints);
			this.spliceWaypoints.apply(this, args);
			return this;
		},

		spliceWaypoints: function() {
			var args = [arguments[0], arguments[1]],
			    i;

			for (i = 2; i < arguments.length; i++) {
				args.push(arguments[i] && arguments[i].hasOwnProperty('latLng') ? arguments[i] : new Waypoint(arguments[i]));
			}

			[].splice.apply(this._waypoints, args);

			// Make sure there's always at least two waypoints
			while (this._waypoints.length < 2) {
				this.spliceWaypoints(this._waypoints.length, 0, null);
			}

			this._updateMarkers();
			this._fireChanged.apply(this, args);
		},

		onAdd: function(map) {
			this._map = map;
			this._updateMarkers();
		},

		onRemove: function() {
			var i;
			this._removeMarkers();

			if (this._newWp) {
				for (i = 0; i < this._newWp.lines.length; i++) {
					this._map.removeLayer(this._newWp.lines[i]);
				}
			}

			delete this._map;
		},

		createGeocoders: function() {
			var container = L.DomUtil.create('div', 'leaflet-routing-geocoders ' + this.options.geocodersClassName),
				waypoints = this._waypoints,
			    addWpBtn,
			    reverseBtn;

			this._geocoderContainer = container;
			this._geocoderElems = [];


			if (this.options.addWaypoints) {
				addWpBtn = L.DomUtil.create('button', 'leaflet-routing-add-waypoint ' + this.options.addButtonClassName, container);
				addWpBtn.setAttribute('type', 'button');
				L.DomEvent.addListener(addWpBtn, 'click', function() {
					this.spliceWaypoints(waypoints.length, 0, null);
				}, this);
			}

			if (this.options.reverseWaypoints) {
				reverseBtn = L.DomUtil.create('button', 'leaflet-routing-reverse-waypoints', container);
				reverseBtn.setAttribute('type', 'button');
				L.DomEvent.addListener(reverseBtn, 'click', function() {
					this._waypoints.reverse();
					this.setWaypoints(this._waypoints);
				}, this);
			}

			this._updateGeocoders();
			this.on('waypointsspliced', this._updateGeocoders);

			return container;
		},

		_createGeocoder: function(i) {
			var geocoder = this.options.createGeocoderElement(this._waypoints[i], i, this._waypoints.length, this.options);
			geocoder
			.on('delete', function() {
				if (i > 0 || this._waypoints.length > 2) {
					this.spliceWaypoints(i, 1);
				} else {
					this.spliceWaypoints(i, 1, new Waypoint());
				}
			}, this)
			.on('geocoded', function(e) {
				this._updateMarkers();
				this._fireChanged();
				this._focusGeocoder(i + 1);
				this.fire('waypointgeocoded', {
					waypointIndex: i,
					waypoint: e.waypoint
				});
			}, this)
			.on('reversegeocoded', function(e) {
				this.fire('waypointgeocoded', {
					waypointIndex: i,
					waypoint: e.waypoint
				});
			}, this);

			return geocoder;
		},

		_updateGeocoders: function() {
			var elems = [],
				i,
			    geocoderElem;

			for (i = 0; i < this._geocoderElems.length; i++) {
				this._geocoderContainer.removeChild(this._geocoderElems[i].getContainer());
			}

			for (i = this._waypoints.length - 1; i >= 0; i--) {
				geocoderElem = this._createGeocoder(i);
				this._geocoderContainer.insertBefore(geocoderElem.getContainer(), this._geocoderContainer.firstChild);
				elems.push(geocoderElem);
			}

			this._geocoderElems = elems.reverse();
		},

		_removeMarkers: function() {
			var i;
			if (this._markers) {
				for (i = 0; i < this._markers.length; i++) {
					if (this._markers[i]) {
						this._map.removeLayer(this._markers[i]);
					}
				}
			}
			this._markers = [];
		},

		_updateMarkers: function() {
			var i,
			    m;

			if (!this._map) {
				return;
			}

			this._removeMarkers();

			for (i = 0; i < this._waypoints.length; i++) {
				if (this._waypoints[i].latLng) {
					m = this.options.createMarker(i, this._waypoints[i], this._waypoints.length);
					if (m) {
						m.addTo(this._map);
						if (this.options.draggableWaypoints) {
							this._hookWaypointEvents(m, i);
						}
					}
				} else {
					m = null;
				}
				this._markers.push(m);
			}
		},

		_fireChanged: function() {
			this.fire('waypointschanged', {waypoints: this.getWaypoints()});

			if (arguments.length >= 2) {
				this.fire('waypointsspliced', {
					index: Array.prototype.shift.call(arguments),
					nRemoved: Array.prototype.shift.call(arguments),
					added: arguments
				});
			}
		},

		_hookWaypointEvents: function(m, i, trackMouseMove) {
			var eventLatLng = function(e) {
					return trackMouseMove ? e.latlng : e.target.getLatLng();
				},
				dragStart = L.bind(function(e) {
					this.fire('waypointdragstart', {index: i, latlng: eventLatLng(e)});
					//this._fireChanged();//μπορω να το βαλώ και εδω
				}, this),
				drag = L.bind(function(e) {
					this._waypoints[i].latLng = eventLatLng(e);
					this.fire('waypointdrag', {index: i, latlng: eventLatLng(e)});
					//this._fireChanged();//μπορω να το βαλώ και εδω
					
					
				}, this),
				dragEnd = L.bind(function(e) {
					this._waypoints[i].latLng = eventLatLng(e);
					this._waypoints[i].name = '';
					if (this._geocoderElems) {
						this._geocoderElems[i].update(true);
					}
					this.fire('waypointdragend', {index: i, latlng: eventLatLng(e)});
					this._fireChanged();
				}, this),
				mouseMove,
				mouseUp;

			if (trackMouseMove) {
				mouseMove = L.bind(function(e) {
					this._markers[i].setLatLng(e.latlng);
					drag(e);
				}, this);
				mouseUp = L.bind(function(e) {
					this._map.dragging.enable();
					this._map.off('mouseup', mouseUp);
					this._map.off('mousemove', mouseMove);
					dragEnd(e);
				}, this);
				this._map.dragging.disable();
				this._map.on('mousemove', mouseMove);
				this._map.on('mouseup', mouseUp);
				dragStart({latlng: this._waypoints[i].latLng});
			} else {
				m.on('dragstart', dragStart);
				m.on('drag', drag);
				m.on('dragend', dragEnd);
			}
		},

		dragNewWaypoint: function(e) {
			var newWpIndex = e.afterIndex + 1;
			if (this.options.routeWhileDragging) {
				this.spliceWaypoints(newWpIndex, 0, e.latlng);
				this._hookWaypointEvents(this._markers[newWpIndex], newWpIndex, true);
			} else {
				this._dragNewWaypoint(newWpIndex, e.latlng);
			}
		},

		_dragNewWaypoint: function(newWpIndex, initialLatLng) {
			var wp = new Waypoint(initialLatLng),
				prevWp = this._waypoints[newWpIndex - 1],
				nextWp = this._waypoints[newWpIndex],
				marker = this.options.createMarker(newWpIndex, wp, this._waypoints.length + 1),
				lines = [],
				draggingEnabled = this._map.dragging.enabled(),
				mouseMove = L.bind(function(e) {
					var i,
						latLngs;
					if (marker) {
						marker.setLatLng(e.latlng);
					}
					for (i = 0; i < lines.length; i++) {
						latLngs = lines[i].getLatLngs();
						latLngs.splice(1, 1, e.latlng);
						lines[i].setLatLngs(latLngs);
					}

					L.DomEvent.stop(e);
				}, this),
				mouseUp = L.bind(function(e) {
					var i;
					if (marker) {
						this._map.removeLayer(marker);
					}
					for (i = 0; i < lines.length; i++) {
						this._map.removeLayer(lines[i]);
					}
					this._map.off('mousemove', mouseMove);
					this._map.off('mouseup', mouseUp);
					this.spliceWaypoints(newWpIndex, 0, e.latlng);
					if (draggingEnabled) {
						this._map.dragging.enable();
					}

					L.DomEvent.stop(e);
				}, this),
				i;

			if (marker) {
				marker.addTo(this._map);
			}

			for (i = 0; i < this.options.dragStyles.length; i++) {
				lines.push(L.polyline([prevWp.latLng, initialLatLng, nextWp.latLng],
					this.options.dragStyles[i]).addTo(this._map));
			}

			if (draggingEnabled) {
				this._map.dragging.disable();
			}

			this._map.on('mousemove', mouseMove);
			this._map.on('mouseup', mouseUp);
		},

		_focusGeocoder: function(i) {
			if (this._geocoderElems[i]) {
				this._geocoderElems[i].focus();
			} else {
				document.activeElement.blur();
			}
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{"./geocoder-element":52,"./waypoint":61}],61:[function(_dereq_,module,exports){
(function (global){
(function() {
	//'use strict';

	var L = (typeof window !== "undefined" ? window['L'] : typeof global !== "undefined" ? global['L'] : null);

	module.exports = L.Class.extend({
		options: {
			allowUTurn: false,
		},
		initialize: function(latLng, name, options) {
			L.Util.setOptions(this, options);
			this.latLng = L.latLng(latLng);
			this.name = name;
		}
	});
})();

}).call(this,typeof global !== "undefined" ? global : typeof self !== "undefined" ? self : typeof window !== "undefined" ? window : {})
},{}]},{},[53]);


//end routing





























function lineDistance(geojson, units) {
    // Input Validation
    if (!geojson) throw new Error('geojson is required');

    // Calculate distance from 2-vertex line segements
    return segmentReduce(geojson, function (previousValue, segment) {
        var coords = segment.geometry.coordinates;
        return previousValue + distance(coords[0], coords[1], units);
    }, 0);
};









































































function helperfeature(geometry, properties, bbox, id) {
    if (geometry === undefined) throw new Error('geometry is required');
    if (properties && properties.constructor !== Object) throw new Error('properties must be an Object');
   // if (bbox && bbox.length !== 4) throw new Error('bbox must be an Array of 4 numbers');
    if (id && ['string', 'number'].indexOf(typeof id) === -1) throw new Error('id must be a number or a string');

    var feat = {type: 'Feature'};
    if (id) feat.id = id;
    if (bbox) feat.bbox = bbox;
    feat.properties = properties || {};
    feat.geometry = geometry;
    return feat;
}





function geometry(type, coordinates, bbox) {
    // Validation
    if (!type) throw new Error('type is required');
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');
    if (bbox && bbox.length !== 4) throw new Error('bbox must be an Array of 4 numbers');

    var geom;
    switch (type) {
    case 'Point': geom = point(coordinates).geometry; break;
    case 'LineString': geom = lineString(coordinates).geometry; break;
    case 'Polygon': geom = polygon(coordinates).geometry; break;
    case 'MultiPoint': geom = multiPoint(coordinates).geometry; break;
    case 'MultiLineString': geom = multiLineString(coordinates).geometry; break;
    case 'MultiPolygon': geom = multiPolygon(coordinates).geometry; break;
    default: throw new Error(type + ' is invalid');
    }
    if (bbox) geom.bbox = bbox;
    return geom;
}













var nativeIsArray = Array.isArray
var toString = Object.prototype.toString;



function isArray(obj) {
    return toString.call(obj) === "[object Array]"
}



var highwaySpeeds = {
    motorway: 110,
    trunk: 90,
    primary: 80,
    secondary: 70,
    tertiary: 50,
    unclassified: 50,
    road: 50,
    residential: 30,  //residential ειναι αυτο που μετραει
    service: 30,
    living_street: 20,
	service: 40,
	asphalt: 50,   //δεν ειμαι σιγουρος για την ταχύτητα
	track:40,	 //δεν ειμαι σιγουρος για την ταχύτητα
	construction: 30,   //δεν ειμαι σιγουρος για την ταχύτητα
	footway: 10,     //δεν ειμαι σιγουρος για την ταχύτητα
	concreteRoad :50,  //δεν ειμαι σιγουρος για την ταχύτητα
	way:40,       //δεν ειμαι σιγουρος για την ταχύτητα
	steps:10,     //δεν ειμαι σιγουρος για την ταχύτητα
	cycleway: 20,     //δεν ειμαι σιγουρος για την ταχύτητα
	path:8,        //δεν ειμαι σιγουρος για την ταχύτητα
    service: 5, //δεν ειμαι σιγουρος για την ταχύτητα
	parking_aisle: 5 , //δεν ειμαι σιγουρος για την ταχύτητα
	living_street :10 , //δεν ειμαι σιγουρος για την ταχύτητα
	
}

var unknowns = {};

/*
function weightFn(a, b, props) {
var dx = a[0] - b[0];
var dy = a[1] - b[1];
return Math.sqrt(dx * dx + dy * dy);
}
*/



function weightFn(a, b, props) {
	//console.log(props.tags.highway,"props.props.tags.highway");

	if (props && props.tags && props.tags.hasOwnProperty("highway")) {
		
	//na γραψω if feature properties contains highwaySpeeds else να αλλαξω το weightFn
    var d = distance(point(a), point(b)) * 1000,
        factor = 0.9,
        type = props.highway || props.tags.highway || props.index.highway || props.index  || props.description || props.type || props.name || props.description.Type ||  props.name.Type ,// εγώ το  άλλαξα
		//type = props.highway || props.tags.highway,
        forwardSpeed,
        backwardSpeed;
 //console.log("d",d);

 
    if (props.maxspeed) {
        forwardSpeed = backwardSpeed = Number(props.maxspeed);
    } 
	
	
	else {
        var linkIndex = type.indexOf('_link');
        if (linkIndex >= 0) {
            type = type.substring(0, linkIndex);
            factor *= 0.7;
        }

        forwardSpeed = backwardSpeed = highwaySpeeds[type] * factor;
        if (!forwardSpeed) {
            unknowns[type] = true;
        }
    }




    if (props.oneway && props.oneway !== 'no' || props.junction && props.junction === 'roundabout') {
        backwardSpeed = null;
    }

    return {
        forward: forwardSpeed && (d / (forwardSpeed / 3.6)),
        backward: backwardSpeed && (d / (backwardSpeed / 3.6)),
    };
	
	
	///else να αλλαξω το weightFn
	
	}
	else
	{
var dx = a[0] - b[0];
var dy = a[1] - b[1];
return Math.sqrt(dx * dx + dy * dy);
		
	}
	
}





var turffeaturecollection = function(e){return{type:"FeatureCollection",features:e}};

var nearest = function(targetPoint, points){
  var nearestPoint;
  var count = 0;
  var dist = Infinity;
  points.features.forEach(function(pt){
    if(!nearestPoint){
      nearestPoint = pt;
      var dist = distance(targetPoint, pt, 'miles');
      nearestPoint.properties.distance = dist;
    }
    else{
      var dist = distance(targetPoint, pt, 'miles');
      if(dist < nearestPoint.properties.distance){
        nearestPoint = pt;
        nearestPoint.properties.distance = dist;
      }
    }
  });
  delete nearestPoint.properties.distance;
  return nearestPoint;
}

 


//PathFinder.extend =

 var rooter = L.Class.extend({
    initialize: function(geojson) {
		//console.log("geojson",geojson);
        this._pathFinder = new PathFinder(geojson, {
            precision: 1e-9,
            weightFn: weightFn
			
        });
        var vertices = this._pathFinder._graph.vertices;

        this._points = turffeaturecollection(Object.keys(vertices)
            .filter(function(nodeName) {
                return Object.keys(vertices[nodeName]).length;
            })
            .map(function(nodeName) {
                var vertice = this._pathFinder._graph.sourceVertices[nodeName];
                return point(vertice);
            }.bind(this)));
        // console.log(JSON.stringify(unknowns, null, 2));
    },

    route: function(waypoints, cb, context) {
        var actualWaypoints = waypoints.map(function(wp) {
                return nearest(utiltoPoint(wp), this._points);
            }.bind(this)),
            legs = actualWaypoints.map(function(wp, i, wps) {
            if (i > 0) {
                return this._pathFinder.findPath(wps[i - 1], wp);
            }

            return [];
        }.bind(this)).slice(1);

        if (legs.some(function(l) { return !l; })) {
            return cb.call(context, {
                status: 1,
                message: 'Can\'t find route.'
            });
        }

        var totalTime = legs.reduce(function(sum, l) { return sum + l.weight; }, 0);
        var totalDistance = legs.reduce(function(sum, l) { 
            var legDistance = l.path.reduce(function(d, c, i, cs) {
                if (i > 0) {
                    return d + distance(point(cs[i - 1]), point(c)) * 1000;
                }
                return d;
            }, 0);
            return sum + legDistance;
        }, 0);

        cb.call(context, null, [{
            name: '',
            waypoints: actualWaypoints.map(function(p) { return { latLng: utiltoLatLng(p) }; }),
            inputWaypoints: waypoints,
            summary: {
                totalDistance: totalDistance,
                totalTime: totalTime
            },
            coordinates: Array.prototype.concat.apply([], legs.map(function(l) { return l.path.map(utiltoLatLng); })),
            instructions: []
        }]);
    }
	
});





    function utiltoPoint (wp) {
        var c = wp.latLng;
        return {
            type: 'Feature',
            geometry: {
                type: 'Point',
                coordinates: [c.lng, c.lat]
            }
        };
    }

    function utiltoLatLng (p) {
        return L.latLng(p[1], p[0]);
    }
