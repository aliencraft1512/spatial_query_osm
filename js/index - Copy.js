

        $("#nearestfunction").show();
            $("#voronoifunction").show();
            $("#pathnetworkfunction").show();
$("#prepareheatmapdata").show();

let selLayer;
let SelectPoints;
let identifyqueryswitch = false;
let esriqueryswitch = false;
let overpassqueryswitch = false;
let southWest; //
let quadtree;
let hash;
let oms;
let omssecond;
let maxzoom = [22];
let sidebar;
let loading;
let individualinfoarray = [];
let selPts = [];
let theidentifymarker;
let boundsrectangle;
let esriquerylayer;
let lat;
let lon;
let xy;
let identify_pts = L.geoJson(null);
let pinned_layer_ = L.geoJson(null);
let $osmrequest;
let within;
let layersarray = [];
let randomlayer = L.featureGroup();
let pointfordistance = null;
let distancepolyline = null;
let differencemarker = null;
let turffeatures = null;
let differencemarkerarray = [];
let voronoifeaturegroup = new L.FeatureGroup();
let voronoipointLayer = new L.geoJSON(null, {
    onEachFeature: function(feature, layer) {
        layer.on('mouseover', function(e) {

            layer.bindTooltip("" + feature.properties.index + "", {
                direction: 'auto',
                sticky: true,
            }).openTooltip();
        })
    },
    pointToLayer: function(feature, latlng) {
        return L.circleMarker(latlng, {
            //renderer: myRenderer,
            radius: 6,
            fillColor: "yellow",
            color: "#000",
            weight: 1,
            opacity: 0,
            fillOpacity: 0,
            boostType: boostType,
            boostScale: boostScale,
            boostExp: boostExp
        })

    }

});


let ptsWithin;
let voronoigeojsoned;
let voronoifeatures;
let voronoiPolygons;
let voronoipointTempLayer;
let topojsonnamearray = [];

//////////////////////////network path geojson oath finder analysis/////////////////
let cleanednetworkarray = [];
let bbox;
let bounds;
let router;
let routingcontrol;
let totalDistance;
let infoContainer;
let networkLayer;
let c1;
let c2;
let node;
let graph;
let nodeNames;
let totalNodes;
let totalEdges;
let controlpoint1;
let controlpoint2;
let network;
//////////////////////////network path geojson oath finder analysis/////////////////



var overlayLayers = {
};


let overpassbound = new L.GeoJSON.AJAX(
    "./data/cyprusbound/overpassbound.geojson", {
        style: {
            opacity: 0,
            fillOpacity: 0
        }
    }
);




let spinopts = {
    lines: 18,
    // The number of lines to draw
    length: 58,
    // The length of each line
    width: 15,
    // The line thickness
    radius: 58,
    // The radius of the inner circle
    scale: 1.05,
    // Scales overall size of the spinner
    corners: 3,
    // Corner roundness (0..1)
    color: "yellow",
    // CSS color or array of colors
    fadeColor: "transparent",
    // CSS color or array of colors
    speed: 1.3,
    // Rounds per second
    rotate: 64,
    // The rotation offset
    animation: "spinner-line-shrink",
    // The CSS animation name for the lines
    direction: 1,
    // 1: clockwise, -1: counterclockwise
    zIndex: 2e9,
    // The z-index (defaults to 2000000000)
    className: "spinner",
    // The CSS class to assign to the spinner
    top: "50%",
    // Top position relative to parent
    left: "50%",
    // Left position relative to parent
    shadow: "0 0 1px transparent",
    // Box-shadow for the lines
    position: "absolute" // Element positioning
};



function hideSpinner() {
    // Hide the spinner
    map.spin(false);
}

function showSpinner() {
    map.spin(true, spinopts);
}





function getvaluefromslider() {
    let val = document.getElementById('Rectangle_value').value;
    return val;
}

function updateTextInput(val) {
    $("#Rectangle_valuelabel").text((val / 1000 * val / 1000).toFixed(1) + " Km²");
}


let rectanglevalue = $("#Rectangle_value").val();
$("#Rectangle_valuelabel").text((rectanglevalue / 1000 * rectanglevalue / 1000).toFixed(1) + " Km²");




let boostType = "circle"; //'balloon';
//'ball';

let boostScale = 1;
let boostExp = 0.125;



let map = new L.Map("map", {

    preferCanvas: true,
    center: [35.1590,33.3637],
    zoom: 16,
    zoomControl: false


});



let googleterain = L.tileLayer(
    "https://{s}.google.com/vt/lyrs=p&x={x}&y={y}&z={z}", {
        maxZoom: maxzoom,
        //detectRetina: true,
        subdomains: ["mt0", "mt1", "mt2", "mt3"]
    }
).addTo(map);
googleterain.crs = L.CRS.EPSG3857; ///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////



sidebar = L.control.sidebar("sidebar", {
    autoPan: false,
    closeButton: true
});
sidebar.addTo(map);

let sidebar2 = L.control.sidebar("sidebars", {
    autoPan: false,
    position: "topleft",
    closeButton: false
});
sidebar2.addTo(map);
sidebar2.show();





function opentab(evt, tabname) {
    sidebar.hide();
    let i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
        sidebar.hide();
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(tabname).style.display = "inline";
    evt.currentTarget.className += " active";
    document.getElementById("sidebars").style.width = "290px";

}





function info(e) {

    let layer = e.target;
    // console.log(layer);
    let lineinfo = "";
    let table = "<table class='table table-striped table-bordered table-sm' style='margin-bottom: 0px;'>";
    let hiddenProps = ["styleUrl", "styleHash", "FILLOPACITY", "fillopacity", "styleMapHash", "stroke", "strokeopacity", "stroke-opacity",
        "stroke-width", "opacity", "fill", "fill-opacity", "icon", "scale", "coordTimes", "_id_", "areaValue", "AREAVALUE_UOM",
        "AREAVALUE_VOID", "BEGINLIFESPANVERSION", "BEGINLIFESPANVERSION_VOID", "ENDLIFESPANVERSION", "ENDLIFESPANVERSION_VOID", "ID_LOCALID",
        "ID_VERSIONID", "ID_VERSIONID_VOID", "VALIDFROM", "VALIDFROM_VOID", "VALIDTO", "VALIDTO_VOID", "BASPROPUNIT_VOID", "ZONING", "ZONING_VOID",
        "ADMUNIT_VOID", "SHAPE.STAREA()", "SHAPE.STLENGTH()", "IFCID", "OBJECTID", "ENDLIFESPANVERSION_VOID", "ENDLIFESPANVERSION", "BEGINLIFESPANVERSION_VOID",
        "BEGINLIFESPANVERSION", "SHAPE.STArea()", "SHAPE.STLength()", "id_versionId_void", "id_versionId", "endLifespanVersion_void",
        "endLifespanVersion", "beginLifespanVersion_void", "beginLifespanVersion", "basPropUnit_void", "areaValue_uom", "zoning_void", "zoning",
        "validTo_void", "validTo", "validFrom_void", "validFrom", "areaValue_void", "admUnit_void", "id_namespace", "id_localId", "admUnit", "gis_GIS__2", "gis_GIS__1", "gis_GIS_dc"
    ];
    for (let k in layer.feature.properties) {
        if (layer.feature.properties.hasOwnProperty(k) && hiddenProps.indexOf(k) == -1) {
            let v = String(layer.feature.properties[k]);

            table += "<tr><th>" + k.toUpperCase() + "</th><td>" + layer.feature.properties[k] + "</td></tr>";
            lineinfo +=
                "<b>" +
                k.toUpperCase() +
                "</b><br>" +
                layer.feature.properties[k] +
                "<br>" +
                '<hr style="margin:5px 0px;  border-style: dotted; border-color: forestgreen;">';
        }
    }
    table += "</table>";

    mapinfo.innerHTML =
        lineinfo;
    individualinfoarray = [];
    individualinfoarray.push(layer.feature);

    sidebar.show(); //sidebarToggle.state('close-sidebar');
}
/*i changed feaure to feature*/

function updateinfo(feature, layer) {

    feature.properties["_id_"] = L.Util.stamp(layer);

    mapinfo.innerHTML = "";
    oms.addMarker(layer);
    quadtree.add(layer);
    layer.on("click", function(e) {

        }),
        layer.on({
            click: info,
            dblclick: info,
            mousedown: info // mouseover: info
            //load: quadtree.add(layer)
        });
} //////////////////////////////////update info function/////////////////////
///////////////////updateinfo normal markers////////////

function updateinfonormalmakrer(feature, layer) {
    feature.properties["_id_"] = L.Util.stamp(layer);

    mapinfo.innerHTML = "";
    omssecond.addMarker(layer);
    quadtree.add(layer);
    layer.on("click", function(e) {
            map.invalidateSize();
        }),
        layer.on({
            click: info,
            dblclick: info,
            mousedown: info // mouseover: info
            //load: quadtree.add(layer)
        });
} ///////////////////updateinfo normal markers////////////
/////////////////updateinfo sites /////////////////////////




let maploadcallBack = function maploadcallBack() {



    $.getScript("./js/osmtogeojson.js").done(function() {

            console.log('./js/osmtogeojson.js loaded');
        })
        .fail(function(jqxhr, settings, exception) {

            console.log('./js/osmtogeojson.js not loaded');
        });




    $.getScript("./js/path.js").done(function() {

            console.log('./js/path.js loaded');

    



        })
        .fail(function(jqxhr, settings, exception) {

            console.log('./js/path.js not loaded');
        });



    $("#identifeatures").show();
    $("#closebuttonidentify, #identifyaccordion, #countfeatures,  #identifyswitches").show();
    identifycontrol.state("hide-identify");
    choosespatialquery();

};



$("#identifybtn").click(function(){
	
   $("#identifeatures").show();
    $("#closebuttonidentify, #identifyaccordion, #countfeatures,  #identifyswitches").show();
    identifycontrol.state("hide-identify");
    choosespatialquery();
})

$("#spatialanalysisbtn").click(function(){
             identifycontrol.state("show-identify");
                closeidentify()
                map.removeLayer(bridgelayer);
})






setTimeout(maploadcallBack, 1000);




var identifycontrol = new L.easyButton({
    position: "topright",
    states: [{
            stateName: "show-identify",
            icon: '<img src="./images/hand.png" style="width:26px">',
            title: "Spatial Query",
            onClick: function onClick(btn, map) {


                map.addLayer(bridgelayer);

                btn.state("hide-identify");
                $("#identifeatures").show();

                choosespatialquery();


                $("#toolbar #toolbarlayerdiv").empty();



                if ($(".mdl-card--identifytopography").hasClass("is-checked")) {

                    map.invalidateSize();

                    if (map.hasLayer(theidentifymarker)) {
                        map.removeLayer(theidentifymarker);
                        map.removeLayer(theidentifyCircle);
                        distance_from_centerPoint = null;
                    }
                    identify_pts.clearLayers();
                    bridgelayer.clearLayers();
                    selPts = [];

                } else {
                    map.invalidateSize();
                    $(".mdl-card--identifytopography").trigger("click");
                    bridgelayer.clearLayers();
                }

                $("#closebuttonidentify, #identifyaccordion, #countfeatures,  #identifyswitches").show();

                map.invalidateSize();

            }



        },
        {
            stateName: "hide-identify",
            icon: '<img src="./images/identify.png" style="width:26px">',
            title: "Hide identify control",
            onClick: function onClick(btn, map) {
                btn.state("show-identify");
                identifycontrol.state("show-identify");
                closeidentify()

                map.removeLayer(bridgelayer);

            }
        }
    ],
    id: 'identifycontrolid'
}).addTo(map);




$(function() {
    let icons = {
        header: "ui-icon-circle-arrow-e",
        activeHeader: "ui-icon-circle-arrow-s"
    };
    $("#identifyaccordion").accordion({
        active: true,
        animate: 200,
        icons: icons,
        collapsible: true,
        heightStyle: "content"
    });
}); //This section creates a bridge layer using the turf.within() after clicking on the
//choropleth layer




function choosespatialquery() {

 $("#closebuttonidentify, #identifyaccordion, #countfeatures,  #identifyswitches").show();

    $("#querycheckboxmenuaccordion").hide();
    $("#overpassquerylayerSelect")
        .fadeOut(150)
        .fadeIn(150);

    $("#toolbar li").remove();
    selPts = [];
    if (map.hasLayer(bridgelayer)) {
        bridgelayer.clearLayers();
    }
    map.addLayer(overpassbound);
    bridgelayer.clearLayers();
    overpassqueryswitch = true;
    let overpassquerylayer = document.getElementById('overpassquerylayerSelect');
    overpassquerylayer.style.display = 'block';

    $("#querycheckboxmenuaccordion").hide();
   overpassquery = overpassquery();
}


$("#closebuttonidentify").on("click", function() {
    destroyanalysis()
    closeidentify()

})




function pinlayer(layer, name) {
    closeidentify();
    layername = name;
    console.log("layername before", layername);
    layername = new L.geoJson(layer, {
        style: {
            color: "black",
            fillOpacity: 0,
            weight: 4
        },
        onEachFeature: function(feature, layer) {
            feature.properties["_id_"] = L.Util.stamp(layer);
            if (feature.geometry.type === 'Polygon') {

                layer.on('mouseover', function(e) {
                    layer.bindTooltip("" + feature.properties.index + "", {
                        direction: 'auto',
                        sticky: true,
                    }).openTooltip();
                })
            } else if (feature.geometry.type === "Point") {

                layer.on('mouseover', function(e) {
                    layer.bindTooltip("" + feature.properties.index + "", {
                        direction: 'auto',
                        sticky: true,
                    }).openTooltip();
                })
            } else {
                let lineinfo = "";
                let hiddenProps = ["styleUrl"];
                for (let k in feature.properties) {
                    if (layer.feature.properties.hasOwnProperty(k) && hiddenProps.indexOf(k) == -1) {
                        let v = String(feature.properties[k]);
                        lineinfo +=
                            "<b>" +
                            k +
                            "</b><br>" +
                            layer.feature.properties[k] +
                            "<br>" +
                            '<hr style="margin:5px 0px;  border-style: dotted; border-color: forestgreen;">';

                    }
                }

                layer.on('click', function(e) {
                    layer.bindPopup("" + lineinfo + "", {
                        direction: 'auto',
                        sticky: false,
                    }).openPopup();
                    console.log(layer);
                })

            }

        },
        pointToLayer: function(feature, latlng) {
            return L.circleMarker(latlng, {
                //renderer: myRenderer,
                radius: 6,
                fillColor: "yellow",
                color: "blue",
                weight: 1,
                opacity: 0,
                fillOpacity: 0,
                boostType: boostType,
                boostScale: boostScale,
                boostExp: boostExp


            })
        }


    }).addTo(map);

    console.log("layername after", layername);
    overlayLayers[L.Util.stamp(layername)] = layername;
    name = name + layername._leaflet_id;
    addLayertolist(layername, name);
    map.addLayer(layername);
    layersarray.push(name);
    randomlayer.addLayer(layername);


}




overpassquery = function overpassquery() {

    $("#overpassquerylayerSelect")
        .fadeOut(150)
        .fadeIn(150)


    identifyqueryswitch = false;
    esriqueryswitch = false;

    if (overpassqueryswitch == true) {


        overpassbound.on("click", function(e) {
            $("#identifeatures").show();
            $("#identifyaccordion").accordion("activate", false);

            xy = null;
            lat = null;
            lon = null;
            lat = e.latlng.lat;
            lon = e.latlng.lng;
            console.log('overpassbound clicked ');
            map.invalidateSize();
            bridgelayer.clearLayers();
            map.invalidateSize();
            $("#identifeatures").show();
            $("#identifyaccordion").accordion("activate", false);
            document.getElementById("countfeatures").innerHTML = "";
            $("#toolbar li").remove();
            hideSpinner();
            showSpinner();
            bridgelayer.clearLayers();
            $("#toolbar li").remove();
            if (boundsrectangle != undefined) {
                map.removeLayer(boundsrectangle);
                map.invalidateSize();
            }

            let overpassApiUrl;




            xy = [lat, lon]; //center point of circle
            selPts = []; //Reset the array 
            var bounds = L.latLng(xy).toBounds(getvaluefromslider());
            boundsrectangle = L.rectangle(bounds, {
                color: 'blue',
                weight: 2,
                fillOpacity: 0.0,
                interactive: false
            }).addTo(map);

            map.fitBounds(bounds);


            function buildOverpassApiUrl(map, overpassQuery) {
                let bounds = boundsrectangle.getBounds().getSouth() + ',' + boundsrectangle.getBounds().getWest() + ',' + boundsrectangle.getBounds().getNorth() + ',' + boundsrectangle.getBounds().getEast();
                let nodeQuery = 'node[' + overpassQuery + '](' + bounds + ');';
                let wayQuery = 'way[' + overpassQuery + '](' + bounds + ');';
                let relationQuery = 'relation[' + overpassQuery + '](' + bounds + ');';
                let query = '?data=[out:json][timeout:15];(' + nodeQuery + wayQuery + relationQuery + ');out body geom;';
                //var baseUrl = 'https://lz4.overpass-api.de/api/interpreter';
                let baseUrl = 'https://overpass-api.de/api/interpreter';
                let resultUrl = baseUrl + query;
                console.log("resultUrl", resultUrl);
                return resultUrl;

            }



            bridgelayer.addTo(map);
            map.invalidateSize();
            let overpassquerylayer = document.getElementById('overpassquerylayerSelect');
            let queryTextfieldValue = overpassquerylayer.value;
            overpassApiUrl = buildOverpassApiUrl(map, queryTextfieldValue);
            console.log(overpassApiUrl);


            $osmrequest = $.ajax({
                url: overpassApiUrl,
                type: 'GET',
                success: function(osmDataAsJson) {
                    bridgelayer.clearLayers();
                    let resultAsGeojson = osmtogeojson(osmDataAsJson);
                    //console.log(resultAsGeojson);
                    //console.log("resultAsGeojson.features", JSON.stringify(resultAsGeojson.features));
                    bridgelayer.addData(resultAsGeojson.features);
                    map.invalidateSize();
                    selPts.push(resultAsGeojson.features);

                    document.getElementById("closebuttonidentify").style.display = "block";
                    document.getElementById("countfeatures").style.display = "block";
                    document.getElementById("identifyswitches").style.display = "block";
                    let headerCell = document.createElement("th");
                    headerCell.className = "headerCell";
                    let headerCell2 = document.getElementById("countfeatures");
                    headerCell2.innerHTML =
                        '<span style="text-transform: uppercase; font-size:14px;">' +
                        selPts[0].length +
                        " " +
                        "features found" +
                        "</span>";

                    if (selPts[0].length > 0) {
                        $("#countfeatures").append("<span id ='pinlayerosm' > &emsp;<i class='fa fa-map-pin' aria-hidden='true' title='create a new layer from your query'></i></span>");
                        let value = 'pinlayerosm';
                        $("#pinlayerosm").click(function() {
                            pinlayer(resultAsGeojson.features, value);
                            swal("the layer is pinned")
                        });
                    }

                    hideSpinner();
                },
                error: function(data) {
                    hideSpinner();

                    if (data.statusText === "Bad Request") {
                        if (overpassqueryswitch == true) {
                            swal('You made a bad request.Try a structured query');
                            map.invalidateSize();
                        }
                        if ($osmrequest !== null) {
                            $osmrequest = null;
                        }
                    } else if (data.statusText === "Too Many Requests") {
                        if (overpassqueryswitch == true) {
                            swal('Hey!! You made more requests that the server can handle. Please WAIT a few minutes and try again');
                            map.invalidateSize();
                        }

                    } else if (data.statusText === "Gateway Timeout") {
                        if (overpassqueryswitch == true) {
                            swal('The server is probably too busy to handle your request.! Try again later');
                            map.invalidateSize();
                        }

                    }


                } //error: function(data) {

            });
            map.invalidateSize();
        }); //on click

    }
}




let bridgelayer = new L.geoJson(null, {

    pointToLayer: function pointToLayer(feature, latlng) {




        if (feature && feature.properties.hasOwnProperty("tags")) {



            let keys = Object.keys(feature && feature.properties.tags);
            keys.forEach(function(key) {
                let p = feature.properties;
                p.index = feature.properties.name + feature.properties.description;
                index = p.index;

            });

            //  return new L.Marker(latlng, {}); 
            return L.circleMarker(latlng, {
                //renderer: myRenderer,
                radius: 6,
                fillColor: "yellow",
                color: "blue",
                weight: 1,
                opacity: 0,
                fillOpacity: 0,
                boostType: 'circle',
                boostScale: boostScale,
                boostExp: boostExp


            }).bindTooltip(feature.properties.index, {
                permanent: false,
                direction: "left"
            });


            if (feature && feature.geometry.type === "LineString") {

                let keys = Object.keys(feature && feature.properties.tags);
                keys.forEach(function(key) {
                    let p = feature.properties;

                    p.index = feature.properties.tags[key];



                });



            } else if (feature && feature.geometry.type === "Polygon") {
                let keys = Object.keys(feature && feature.properties.tags);
                keys.forEach(function(key) {
                    let p = feature.properties;
                    p.index = feature.properties.tags[key];

                });

            }




        } else if (Object.values(feature.properties).includes("null")) {
            let toponymicon = L.icon({
                iconUrl: "./images/toponym.png",
                //shadowUrl: 'leaf-shadow.png',
                iconSize: [10, 10],
                // size of the icon
                className: "toponymicon"
            });
            return new L.Marker(latlng, {
                icon: toponymicon
            });
        } else if (Object.values(feature.properties).includes("")) {

            return L.circleMarker(latlng, {
                color: "#fff",
                stroke: false,
                weight: 0,
                fillOpacity: 0.4,
                fillColor: "yellow",
                radius: 7,
                boostType: "balloon",
                boostScale: boostScale,
                boostExp: boostExp,
                bubblingMouseEvents: false
            })

        } else {


            return L.circleMarker(latlng, {
                color: "#fff",
                stroke: false,
                weight: 0,
                fillOpacity: 0.4,
                fillColor: "yellow",
                radius: 7,
                boostType: "balloon",
                boostScale: boostScale,
                boostExp: boostExp,
                bubblingMouseEvents: false
            })



        }

    },



    onEachFeature: function onEachFeature(feature, layer) {

  let popupContent = "";
        if (layer.feature && layer.feature.properties.hasOwnProperty("label")) {
            let p = feature.properties;
            p.index = feature.properties.label;


        } else if (layer.feature && layer.feature.properties.hasOwnProperty("tags")) {

            let keys = Object.keys(layer.feature && layer.feature.properties.tags);
            keys.forEach(function(key) {
                let p = feature.properties;
                p.index = feature.properties.name + feature.properties.description;

                index = p.index;



            });



            if (layer.feature && layer.feature.geometry.type === 'LineString') {


                let keys = Object.keys(layer.feature && layer.feature.properties.tags);
                keys.forEach(function(key) {
                    let p = feature.properties;

                    if (layer.feature && layer.feature.properties.tags.hasOwnProperty("name:en:")) {

                        p.index = "Type:  " + feature.properties.tags[key] + "Road Name: " + feature.properties.tags["name:en:"];

                        layer.setStyle({
                            color: 'red',
                            weight: 3,
                        });
						
                    } else if (layer.feature && layer.feature.properties.tags.hasOwnProperty("name")) {

                        p.index = "Type:  " + feature.properties.tags[key] + "Road Name: " + feature.properties.tags.name;

                        layer.setStyle({
                            color: 'red',
                            weight: 3,
                        });
						
                    } 
					else
                    {
                        p.index = "Type:  " + feature.properties.tags[key];
                        // p.tags = feature.properties.tags[key];
                        layer.setStyle({
                            color: 'red',
                            weight: 3,
                        });


                    }


                    layer.setStyle({
                        color: 'red',
                        weight: 6,

                    });


                });




            } else if (layer.feature && layer.feature.geometry.type === "Polygon") {


                let keys = Object.keys(layer.feature && layer.feature.properties.tags);
                keys.forEach(function(key) {
                    let p = feature.properties;
                    p.index = feature.properties.tags[key];
					 popupContent = "";
                    popupContent = popupContent + "<dt>" + key + "</dt><dd>" + feature.properties.tags[key] + "</dd>";


                });


                layer.setStyle({
                    color: 'orange',
                    weight: 6,

                });


            }




        } else {
			
			      let p = feature.properties;
        p.name = feature.properties.index;
        p.description = feature.properties.index;




        layer.feature.id = feature.properties.index;
			
		}




  
        //layer._leaflet_id = feature.properties.index;


        if (layer.feature && layer.feature.properties.hasOwnProperty("tags")) {
            p = layer.feature.properties;
            p.index = layer.feature.properties.name + layer.feature.properties.description;
            index = p.index;
          
            popupContent = popupContent + "<dt>@id</dt><dd>" + feature.properties.type + "/" + feature.properties.id + "</dd>";
            let keys = Object.keys(feature.properties.tags);
            keys.forEach(function(key) {
                popupContent = popupContent + "<dt>" + key + "</dt><dd>" + feature.properties.tags[key] + "</dd>";
            });
            popupContent = popupContent + "</dl>"
            layer.bindPopup(popupContent).bindTooltip(popupContent, {
                direction: 'auto',
                sticky: true,
            });
        } else

        {
			
			p = layer.feature.properties;
            p.index = layer.feature.properties.name + layer.feature.properties.description;
            index = p.index;
			
            popupContent = "" + feature.properties.index;
            let identifypopup = layer.bindPopup(popupContent);

        }

        layer.on({
            click: info
        }); 
 
        layer.on("dblclick", function() {
            //layer.bindPopup(popupContent);
            identifypopup.openPopup();
        }); //quadtree.add(layer);


    }
}).addTo(map);




///////////turf analysis functions /////////////////


function createradiomenu() {
 $("#leftpanel").show();
    destroyanalysis();
    closeidentify();
    $("#closebuttonidentify").show();
    $("#identifeatures").show();
    $("#countfeatures").show();
    $("#countfeatures").html("Choose Layer to Analyze");

    layersarray.forEach(function(item, index) {
        console.log(item);

        let values = item;
        let confirmbutton = document.createElement("button");
        confirmbutton.textContent = "Confirm";
        let queryel = document.createElement("INPUT");
        queryel.setAttribute("type", "radio");
        queryel.setAttribute('style', 'margin-left: 8px!important;');
        //queryel.type = "checkbox";
        queryel.name = "layerid";
        queryel.id = "layerid";
        queryel.value = values;
        queryel.textContent = "";
        let hr = document.createElement('hr');
        let label = document.createElement('span');
        label.htmlFor = "layerid";
        label.id = "labellayerid";
        label.appendChild(document.createTextNode(values));
        let overall = document.createElement("div");
		

		
		
        overall.id = "toolbarlayerdiv";
        overall.appendChild(hr);
        overall.appendChild(label);
        overall.appendChild(queryel);
        let selectdiv = document.getElementById('toolbar');
        selectdiv.appendChild(overall);
		  $("#toolbarlayerdiv")
        .fadeOut(150)
        .fadeIn(150)
		 .fadeOut(150)
        .fadeIn(150)
		 .fadeOut(150)
        .fadeIn(150);

    })


}




function choosepointlayer() {

	  $(".closetabs").trigger("click");
    createradiomenu()
    let value;
	 $("#leftpanel").show();
    $("#leftpanel input").change(function() {
        //value= $("input:radio:checked").val();

        value = $('input:checked', '#leftpanel').val()

        swal({
            title: "!!!!",
            text: "Confirm",
            icon: "info",
            buttons: ["Confirm", "cancel"],
            dangerMode: true
        }).then(function(isConfirm) {
            if (isConfirm) {
                destroyanalysis()

            } else {
                console.log("value", value);
                createturffeaturecollection(value);
                value = null;
                $("#toolbar #toolbarlayerdiv").empty();
                //$("#leftpanel").hide();
                $("#countfeatures").html("Results");



            }
        }); //swal({title: "!!!!",


    });


}




function createturffeaturecollection(sentlayer) {
    differencemarkerarray = [];
    let layer;
    let analysedlayer;
    if (sentlayer.includes('pinlayer')) {
        let layerid = sentlayer.replace(/^\D+/g, '');
        console.log(layerid);
        layer = randomlayer.getLayer(layerid);
        analysedlayer = layer.toGeoJSON();
        console.log("toGeoJSON()", layer);

    } else {
        layer = window[sentlayer];
        console.log("window[sentlayer]", layer);
        layer = layer.toGeoJSON();
        console.log("layer.toGeoJSON()", layer);
        analysedlayer = layer;
        layer = null;
    }

    try {

        let features = turf.featureCollection(analysedlayer.features);
        console.log("features", features);
        differencemarker = new L.Marker(null, {
            draggable: true
        });
        differencemarker.setLatLng(map.getCenter())
        differencemarker.addTo(map).bindPopup("drag marker to find its proximity to your selected features").openPopup();
        let differencemarkerarray = [differencemarker._latlng.lng, differencemarker._latlng.lat];
        let turfpointfinal = turf.point([differencemarker._latlng.lng, differencemarker._latlng.lat], {
            "marker-color": "#0F0"
        });
        let nearest = turf.nearestPoint(turfpointfinal, features);
        nearestfunction(features, differencemarkerarray)
        differencemarker.on("dragend", function(event) {
            //differencemarker.bindPopup("drag marker around to see the proximity of features to it").openPopup()
            differencemarkerarray = [differencemarker._latlng.lng, differencemarker._latlng.lat];
            nearestfunction(features, differencemarkerarray)
        });

    } catch (error)

    {

        if (error.toString().indexOf('coord must be GeoJSON Point or an Array of numbers')) {
            swal("Wrong layer.This layer does not consist of Points");
            destroyanalysis()
        }
        layer = null;
        analysedlayer = null;
        destroyanalysis()
    }
    layer = null;
    analysedlayer = null;
}




function nearestfunction(features, turfpoint) {

    turffeatures = null;
    $("#toolbar ul").html("");
    if (pointfordistance != undefined) {
        map.removeLayer(pointfordistance);
    }
    if (distancepolyline != undefined) {
        map.removeLayer(distancepolyline);
    }
    let redIcon = new L.Icon({
        iconUrl: 'https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-2x-green.png',
        shadowUrl: 'https://cdnjs.cloudflare.com/ajax/libs/leaflet/0.7.7/images/marker-shadow.png',
        iconSize: [25, 41],
        iconAnchor: [12, 41],
        popupAnchor: [1, -34],
        shadowSize: [41, 41]
    });

    pointfordistance = L.geoJSON(null, {
        pointToLayer: function(feature, latlng) {
            return L.marker(latlng, {
                icon: redIcon
            });
        }
    }).addTo(map);



    turffeatures = features;
    let turfpointfinal = turf.point(turfpoint, {
        "marker-color": "#0F0"
    });
    let nearest = turf.nearestPoint(turfpointfinal, turffeatures);
    pointfordistance.addData(nearest);


    distancepolyline = L.polyline([
        [turfpoint[1], turfpoint[0]],
        [nearest.geometry.coordinates[1], nearest.geometry.coordinates[0]]
    ], {
        color: 'yellow',
        weight: 3,
        opacity: 1,
        dashArray: "5, 3"
    }).addTo(map);
    let from = turf.point([turfpoint[1], turfpoint[0]]);
    let to = turf.point([nearest.geometry.coordinates[1], nearest.geometry.coordinates[0]]);
    let options = {
        units: 'kilometers'
    };
    let distance = turf.distance(from, to, options);
    //var distance = turfpoint.distanceTo(turf.getCoord(nearest))/1000;
    //$("#toolbar ul").html(" "+nearest.properties.index+" Position  is "+distance.toFixed(3)+" km away ");
    let $listItem = $("<li>").html(" " + nearest.properties.index + " Position  is " + distance.toFixed(3) + " km away ").appendTo("#toolbar ul");
    distancepolyline.bindTooltip("" + distance.toFixed(3) + " km away ", {
        direction: 'auto',
        sticky: true,
    })
    pointfordistance.bindPopup(" " + nearest.properties.index + " Position  is " + distance.toFixed(3) + " km away ", {
        direction: 'auto',
        sticky: true,
    }).openPopup();


}




function destroyanalysis() {
    try {

        differencemarkerarray = [];
        $("#closebuttonidentify").hide();
        if (pointfordistance != undefined) {
            pointfordistance.clearLayers();
            map.removeLayer(pointfordistance);
            pointfordistance = null;
        }
        if (distancepolyline != undefined) {
            map.removeLayer(distancepolyline);
            distancepolyline = null;
        }
        if (differencemarker != undefined) {
            map.removeLayer(differencemarker);
        }

        if (map.hasLayer(voronoifeaturegroup)) {
            voronoifeaturegroup.clearLayers()
            map.removeLayer(voronoifeaturegroup);
        }

        if (map.hasLayer(voronoipointLayer)) {
            voronoipointLayer.clearLayers()
            map.removeLayer(voronoipointLayer);
        }

        ptsWithin = null;
        voronoigeojsoned = null;
        voronoifeatures = null;
        voronoiPolygons = null;
        voronoipointTempLayer = null;
        $("#toolbar #toolbarlayerdiv").empty();
        $("#toolbar li").remove();
        $("#leftpanel").hide();
        $("#countfeatures").html("");
        $("#countfeatures").hide();

        $("#toolbar #toolbarlayerdiv").remove();
        $("#leftpanel input").empty();

        if (routingcontrol != undefined) {

            map.removeControl(routingcontrol);
            routingcontrol = null;
        }
        router = null;

        totalDistance = null;
        infoContainer = null;
        networkLayer = null;
        c1 = null;
        c2 = null;
        node = null;
        graph = null;
        nodeNames = null;
        totalNodes = null;
        totalEdges = null;
        controlpoint1 = null;
        controlpoint2 = null;
        network = null;
        router = null;
    } catch (error)

    {
        console.log("destroy analysis error", error);
    }
}




function voronoi(sentlayer) {
    voronoifeaturegroup.addTo(map);
    voronoipointLayer.addTo(map);

    let layer;
    if (sentlayer.includes('pinlayer')) {
        let layerid = sentlayer.replace(/^\D+/g, '');
        console.log(layerid);
        layer = randomlayer.getLayer(layerid);
        voronoigeojsoned = layer.toGeoJSON();

    } else {
        layer = window[sentlayer];
        voronoigeojsoned = layer.toGeoJSON();

    }




    // minx is the longitude of the southwestern corner

    //miny is the latitude of the southwestern corner

    //maxx is the longitude of the northeastern corner

    //maxy is the latitude of the northeastern corner


    //minX, minY, maxX, maxY

    let layerWEST = layer.getBounds().getWest();
    let layerEAST = layer.getBounds().getEast();
    let layerNORTH = layer.getBounds().getNorth();
    let layerSOUTH = layer.getBounds().getSouth();
    let WEST = map.getBounds().getWest();
    let EAST = map.getBounds().getEast();
    let NORTH = map.getBounds().getNorth();
    let SOUTH = map.getBounds().getSouth();

    layer = null;


    let boundoptions;

    swal({
        title: "!!!!",
        text: "Confirm",
        icon: "info",
        buttons: ["Use current Map bounds", "Use layer extension bounds"],
        dangerMode: true
    }).then(function(isConfirm) {
        if (isConfirm) {

            boundoptions = {
                bbox: [layerWEST, layerSOUTH, layerEAST, layerNORTH]
            };

        } else {

            boundoptions = {
                bbox: [WEST, SOUTH, EAST, NORTH]
            };

        }


        hideSpinner();
        showSpinner();
        setTimeout(function() {

            try {


                voronoifeatures = turf.featureCollection(voronoigeojsoned.features, boundoptions);
                voronoiPolygons = turf.voronoi(voronoifeatures, boundoptions);

            } catch (error) {
                hideSpinner();
                console.log(error);
                if (error.toString().indexOf('must be a Point')) {
                    swal("Wrong layer.This layer does not consist of Points");
                }


            }
            voronoifeaturegroup.clearLayers();
            voronoipointLayer.clearLayers();
            voronoipointTempLayer = new L.geoJSON(voronoigeojsoned.features);

            console.log("voronoiPolygons.length", voronoiPolygons.features.length);

            try {

                turf.featureEach(voronoiPolygons, function(currentFeature, featureIndex) {

                    //console.log(featureIndex);
                    if (currentFeature && currentFeature.geometry && currentFeature.geometry.type === "Polygon") {
                        let flipped = turf.flip(currentFeature);

                        let polygon = L.polygon(flipped.geometry.coordinates, {
                            weight: 4,
                            fillOpacity: 0.07,
                            color: 'black',
                            dashArray: '1',
                            fillColor: '#0000FF',
                            pane: "topoPane"
                        }).addTo(voronoifeaturegroup);

                    }

                });

            } catch (error) {
                console.log(error)
                hideSpinner();
                swal("An error occured");
            }

            ptsWithin = null;
            ptsWithin = turf.within(voronoipointTempLayer.toGeoJSON(), voronoifeaturegroup.toGeoJSON());
            voronoipointLayer.addData(ptsWithin);
            map.fitBounds(voronoifeaturegroup.getBounds());
            ptsWithin = null;
            $("#countfeatures").html('create a new layer');
            $("#countfeatures").append("<span id ='pinlayervoronoi' > &emsp;<i class='fa fa-map-pin' aria-hidden='true' title='create a new layer'></i></span>");
            voronoifeaturegroup.addLayer(voronoipointLayer);
            let value = 'pinlayervoronoi';
            $("#pinlayervoronoi").click(function() {

                map.flyTo([map.getCenter().lat, map.getCenter().lng], 13)
                let temp;
                temp = voronoifeaturegroup.toGeoJSON();
                pinlayer(temp, value);
                swal("the layer is pinned")


                temp = null;
                voronoigeojsoned = null;
                destroyanalysis()
            });



            hideSpinner();
        }, 200);

    }); //swal({title: "!!!!",


}




function choosevoronoilayer() {
	  $(".closetabs").trigger("click");
	  
    createradiomenu()
	 $("#leftpanel").show();
    $("#leftpanel input").change(function() {
        //let value= $("input:radio:checked").val();

        let value = $('input:checked', '#leftpanel').val()


        swal({
            title: "!!!!",
            text: "Confirm",
            icon: "info",
            buttons: ["cancel", "Confirm"],
            dangerMode: true
        }).then(function(isConfirm) {
            if (isConfirm) {
                console.log("value", value);
                voronoi(value)
                $("#toolbar #toolbarlayerdiv").empty();
                //$("#leftpanel").hide();
                $("#countfeatures").html("Results");


            } else {

                destroyanalysis()
            }
        }); //swal({title: "!!!!",

    });


}




function initializenetwork()

{
	$(".closetabs").trigger("click");
    destroyanalysis()
    createradiomenu()
    let value;
	$("#leftpanel").show();
    $("#leftpanel input").change(function() {
        //value = $("input:radio:checked").val();
        value = $('input:checked', '#leftpanel').val();


        swal({
            title: "!!!!",
            text: "Confirm",
            icon: "info",
            buttons: ["cancel", "Confirm"],
            dangerMode: true
        }).then(function(isConfirm) {
            if (isConfirm) {

                $("#toolbar #toolbarlayerdiv").empty();
                //$("#leftpanel").hide();
                $("#countfeatures").html("Results");
                console.log("value", value);
                initializenetworklayer(value)
                value = null;
            } else {
                destroyanalysis()
            }
        });

    });

}




function initializenetworklayer(sentlayer) {
    cleanednetworkarray = [];
    let layer;
    let cleanedlayer;
    if (sentlayer.includes('pinlayer')) {
        let layerid = sentlayer.replace(/^\D+/g, '');
        console.log(layerid);
        layer = randomlayer.getLayer(layerid);
        cleanedlayer = layer;
        layer = null;
    } else {
        layer = window[sentlayer];
        cleanedlayer = layer;
        layer = null;
    }

    cleanedlayer.eachLayer(function(item) {
        //if((item.feature.geometry.type === "Point" )|| (item.feature.geometry.type === "multipolygon")){

        //}
        //else
        //{

        cleanednetworkarray.push(item.feature);

        //}
    });


    /*
     var uncleannetworkarray=[];
     for (var i = 0; i < cleanednetworkarray.length; i++) { 
    if((cleanednetworkarray[i].geometry.type === "Linestring") || (cleanednetworkarray[i].geometry.type === "Polygon" )){
    uncleannetworkarray.push(cleanednetworkarray[i]);
    console.log("uncleannetworkarray",cleanednetworkarray[i]);
    }}
    */


    exportedlayer = {};
    exportedlayer.type = "FeatureCollection";
    exportedlayer.features = {};
    exportedlayer.features = cleanednetworkarray;
    initializenetworkanalysis(exportedlayer);

    cleanedlayer = null;
}




function initializenetworkanalysis(network) {

    console.log("network", network);
    try {

        $("#countfeatures").hide();
        $("#leftpanel").hide();
        $("#closebuttonidentify").hide();
        bbox = extent(network);

        //bounds = L.latLngBounds([bbox[1], bbox[0]], [bbox[3], bbox[2]]);
        // map.fitBounds(bounds);
        //L.rectangle(bounds, {color: 'orange', weight: 1, fillOpacity: 0.03, interactive: false}).addTo(map);
        router = new rooter(network);

        console.log("router", router);


        controlpoint1 = [router._points.features[0].geometry.coordinates[1], router._points.features[0].geometry.coordinates[0]];
        controlpoint2 = [router._points.features[1].geometry.coordinates[1], router._points.features[1].geometry.coordinates[0]];


        routingcontrol = L.Routing.control({
            createMarker: function(i, wp) {
                return L.marker(wp.latLng, {
                    icon: L.icon.glyph({
                        prefix: '',
                        glyph: String.fromCharCode(65 + i)
                    }),
                    draggable: true
                })
            },
            router: router,
            routeWhileDragging: true,
            routeDragInterval: 100
        }).addTo(map);

        routingcontrol.setWaypoints([
            controlpoint1,
            controlpoint2,
        ]);

        totalDistance = network.features.reduce(function(total, feature) {
            if (feature.geometry.type === 'LineString') {
                return total += lineDistance(feature, 'kilometers');
            } else {
                return total;
            }
        }, 0);


        //if(router._graph.compactedVertices){
        graph = router._pathFinder._graph.compactedVertices;
        //}
        nodeNames = Object.keys(graph);
        totalNodes = nodeNames.length;
        totalEdges = nodeNames.reduce(function(total, nodeName) {
            return total + Object.keys(graph[nodeName]).length;
        }, 0);

        infoContainer = document.querySelector('#info-container');

        [
            ['Total Road Length', totalDistance, 'km'],
            ['Network Nodes', totalNodes / 1000, 'k'],
            ['Network Edges', totalEdges / 1000, 'k'] //,
            // ['Coordinates', router._pathFinder._points.features.length / 1000, 'k']
        ].forEach(function(info) {
            var li = L.DomUtil.create('li', '', infoContainer);

            console.log("info", info);
            li.innerHTML = info[0] + ': <strong>' + Math.round(info[1]) + (info[2] ? '&nbsp;' + info[2] : '') + '</strong>';


        });


        var vertices = router._pathFinder._graph.sourceVertices;
        nodeNames.forEach(function(nodeName) {
            node = graph[nodeName];
            Object.keys(node).forEach(function(neighbor) {
                c1 = vertices[nodeName];
                c2 = vertices[neighbor];
                //  L.polyline([[c1[1], c1[0]], [c2[1], c2[0]]], { weight: 1, opacity: 0.4, renderer: renderer, interactive: false })
                //     .addTo(networkLayer)
                //     .bringToBack();
            });
        });

    } catch (error)

    {

        console.log(error)
        if (error.toString() === "TypeError: Cannot read properties of undefined (reading 'compactedVertices')") {
            //swal("!!The layer contains points"); 
        } else if (error.toString() === 'Error: Compacted graph contains no forks (topology has no intersections).') {
            swal("!!Wrong layer.This layer is not a road network but a point Layer");
        }



    } finally {
        totalDistance = null;
        infoContainer = null;
        networkLayer = null;
        c1 = null;
        c2 = null;
        node = null;
        graph = null;
        nodeNames = null;
        totalNodes = null;
        totalEdges = null;
        controlpoint1 = null;
        controlpoint2 = null;
        network = null;
        router = null;

    }

}




function prepareheatmapdata() {

    {
        destroyanalysis()
        createradiomenu()

	$("#countfeatures").append("<span class='heatmapicon' style='background: url(./images/heatmap2.png) no-repeat;float: left;width: 40px;height: 30px; margin-top: -6px;'></span>");
	
	
        $("#leftpanel input").on('change', function(e) {
            //value = $("input:radio:checked").val();
            value = $('input:checked', '#leftpanel').val();


            swal({
                title: "!!!!",
                text: "Confirm",
                icon: "info",
                buttons: ["cancel", "Confirm"],
                dangerMode: true
            }).then(function(isConfirm) {
                if (isConfirm) {

                    $("#toolbar #toolbarlayerdiv").empty();
                    //$("#leftpanel").hide();
					 $("#countfeatures").html("");
                    $("#countfeatures").html("Results");
                    //console.log("value", value);
                    executeheatmap(value)

                } else {
                    destroyanalysis()

                }
            });

        });

    }

}




function executeheatmap(data) {

layername= data;


    if (data.includes('pinlayer')) {
        //let layerid= data.replace( /^\D+/g, '');
        let layerid = data.replace(/(^.*\[|\].*$)/g, '');

        console.log(layerid);
        //data = randomlayer.getLayer(layerid);
        data = overlayLayers[layerid];
        data = data.toGeoJSON();
        console.log("toGeoJSON()", data);

    } else {

        data = window[data];
        data = data.toGeoJSON();

    }


try{
	turf.voronoi(data, null);
}


catch(error)
{

console.log(error);

if  (error.toString() === "Error: Invalid input to points: must be a Point, given Polygon") {
   swal("The layer does not consist of points but of Polygons..!! The data has been converted to the polygons' centroid values!!"); 
   
       
          var centroids = {
                      type: 'FeatureCollection',
                      features: data.features.map(function(feature) {
                          return {
                              type: 'Feature',
                              properties: feature.properties,
                              geometry: turf.centroid(feature).geometry}})};	   
console.log("centroids",centroids);

  data = centroids; 
     } 

}

finally 
{
	


    datafeatures = turf.featureCollection(data.features);
    var heatmapdata = turf.flip(datafeatures);
    console.log(heatmapdata);
    if (heatmapdata.features.length > 0) {
        heatmapdata = heatmapdata.features.map(function(p) {
            return [p.geometry.coordinates[0], p.geometry.coordinates[1]];
        });
        //console.log(heatmapdata);
        heatlayer = L.heatLayer(heatmapdata, {
            minOpacity: 0.5,
            max: 1,
            radius: 25
        }).addTo(map);

var bbox = turf.bbox(datafeatures);//bbox extent in minX, minY, maxX, maxY order		
console.log("bbox",bbox);
   var b0 = bbox[0];
   var b1 = bbox[1];
   var b2 = bbox[2];
   var b3 = bbox[3]; 
		
map.fitBounds([[b1, b0],[b3, b2]]);

        //var radiusheat = L.DomUtil.get('radiusheat');
        //console.log(radiusheat);
        //radiusheat.value= heatlayer.options.radius;

        /*
        console.time('delaunay');
        var delaunayarray=[];
        for (var i = 0; i < heatmapdata.length; i++) {	
        delaunayarray.push(heatmapdata[i][0],heatmapdata[i][1]);
        }
        var delaunay = new Delaunator(delaunayarray);
        console.log("delaunay",delaunay);
        console.timeEnd('delaunay');
*/

$("#heatmapdivpanel").show();
$("#countfeatures").append(' for Layer:<b>'+layername+'</b>');
 if( $('#radiusheat').length )       
{
	
}
else
{
$("#heatmapdivpanel").append('<br><input  type="range" id="radiusheat" style="width:80%" min="10"  max="30"  step="0.01" class="mdl-slider mdl-js-slider" />'); 
$("#heatmapdivpanel .mdl-slider__container").remove();


}

$("#heatmapdivpanel").appendTo("#toolbar");


[].forEach.call(
            document.querySelectorAll("#heatmapdivpanel input"),
            function(input) {
                input["oninput" in input ? "oninput" : "onchange"] = function(e) {
                    updateheatValues();

                };
            }
        );
 //var maxheat = document.getElementById('maxheat');
 } else
 {
 swal("This Layer does not contain any data");
}




}

}







function updateheatValues() {


    //console.log("radiusheatvalue",radiusheat.value);	

    //heatlayer.setOptions({minOpacity:minOpacityvalue.value});
    heatlayer.setOptions({
        radius: radiusheat.value
    });
    //heatlayer.setOptions({max:maxheat.value});
    //heatlayer.setOptions({radius:12});
    map.invalidateSize();

}





















if ($("#map").length) {
    var mapDiv = document.getElementById("map");
    var overlayDiv = document.getElementById("identifeatures");
    var rectanglequery = L.DomUtil.get('Rectangle_size');
    mapDiv.appendChild(overlayDiv);
    var div = L.DomUtil.get('identifeatures');
    if (!L.Browser.touch) {
        L.DomEvent.disableClickPropagation(div);
        L.DomEvent.on(div, 'mousewheel', L.DomEvent.stopPropagation);
        L.DomEvent.disableScrollPropagation(div);
        L.DomEvent.disableClickPropagation(rectanglequery);
        L.DomEvent.on(rectanglequery, 'mousewheel', L.DomEvent.stopPropagation);
        L.DomEvent.disableScrollPropagation(rectanglequery);

    } else {
        L.DomEvent.on(div, 'click', L.DomEvent.stopPropagation);
        L.DomEvent.disableScrollPropagation(div);


        L.DomEvent.on(rectanglequery, 'mouseenter', function(e) {
            map.dragging.disable()
        });
        L.DomEvent.on(rectanglequery, 'mouseleave', function(e) {
            map.dragging.enable();
        });


        $("#coordinateslatlon, #identifeatures").mousedown(function() {
            map.dragging.disable();
        });

        $(document).mouseup(function() {
            map.dragging.enable();

        });



    }


}




function closeidentify() {
    map.removeLayer(bridgelayer);
    $("#identifeatures").hide();
    identifycontrol.state("show-identify");
    $("#querycheckboxmenuaccordion").hide();



    map.removeLayer(overpassbound);
    xy = [];
    lat = [];
    lon = [];
    overpassqueryswitch = false;

    within = null;

    if ($osmrequest !== null) {
        $osmrequest = null;
    }

    if (!selLayer == null) {
        selLayer.clearLayers();
        map.removeLayer(selLayer);
    }
    selLayer = null;

    if (L.Browser.mobile) {
        // identifycontrol.state("show-identify");
        // swal("This function is not available in mobile");
    } else {
        $("#identifyaccordion,#leftpanel,#countfeatures,#closebuttonidentify,#identifyswitches").hide();
        hideSpinner();

        if ($osmrequest != null) {
            $osmrequest.abort();
            $osmrequest = null;

        }

        document.getElementById("countfeatures").innerHTML = "";
        $("#toolbar li").remove();
        identifycontrol.state("show-identify");
        $("#querycheckboxmenuaccordion").hide();
        var headerCell2 = document.getElementById("countfeatures");
        headerCell2.innerHTML = "";
        headerCell2.innerHTML += "click map to uncover features";
        bridgelayer.clearLayers();
        identify_pts.clearLayers();

    }

    if (boundsrectangle != undefined) {
        map.removeLayer(boundsrectangle);
    }


}




let closesidebarsss = document.querySelector(".closetabs");
closesidebarsss.onclick = function() {
    document.getElementById("spatialanalysis").style.display = "none";
	//$("#identifybtn").trigger("click");
    sidebar.hide();
};





let layerControl = new L.control.layers(null, overlayLayers, {
    position: "topright",
    collapsed: true,
    autoZIndex: true
});



function addLayertolist(layer, name, identifier) {


    layerControl.addOverlay(layer, `
    <span class="layer-name" id="${L.Util.stamp(layer)}">
      ${name}
    </span>
    <span class="layer-buttons">
      <span style="display: ${layer instanceof L.GeoJSON ? 'none' : 'unset'}">
        <a class="layer-btn" href="#" title="Change opacity" onclick="changeOpacity(${L.Util.stamp(layer)}); return false;"><i class="fas fa-adjust"></i></a>
      </span>
	  
	
		  

      <a class="layer-btn" href="#" title="Zoom to layer" onclick="zoomToLayer(${L.Util.stamp(layer)}); return false;"><i class="fas fa-expand-arrows-alt"></i></a>
      <a class="layer-btn" href="#" title="Remove layer" onclick="removeLayer(${L.Util.stamp(layer)}, '${name}'); return false;"><i class="fas fa-trash" style="color: red"></i></a>
    </span>
    <div style="clear: both;"></div>
  `);

    if (layer instanceof L.GeoJSON) {
       // layersarray.push('pinlayer: ' + name + '[' + L.Util.stamp(layer) + ']');
    } else {

    }
    map.invalidateSize();
}
