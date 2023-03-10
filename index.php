




<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
 "http://www.w3.org/TR/html4/strict.dtd">
 

<html lang="en">

<head>
						  
<title>CMPmap</title>



<meta charset="UTF-8">
<!--meta charset="utf-8" /-->
    <meta name="format-detection" content="telephone=no" />
    <meta name="msapplication-tap-highlight" content="no" />

	  <meta http-equiv="Content-Security-Policy" content="upgrade-insecure-requests"  />
    <meta name="viewport" content="user-scalable=no, initial-scale=1, maximum-scale=1, minimum-scale=1, width=device-width" />
    <!-- This is a wide open CSP declaration. To lock this down for production, see below. -->
    <!--meta http-equiv="Content-Security-Policy" content="default-src * 'unsafe-inline'; style-src 'self' 'unsafe-inline'; media-src *" /-->
 
  <meta name="mobile-web-app-capable" content="yes">
    <meta name="apple-mobile-web-app-capable" content="yes">
<meta name="theme-color" content="#263238"/>

   <meta name="apple-mobile-web-app-status-bar-style" content="black">
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />






<style>



		
html, body {
    height: 100%;
    width: 100vw;
	margin:0!important;
	padding:0!important;
}
		

 #map {
        /*  position: absolute;
       height: 100vh;*/
	   position: absolute;
         width: 100vw;
		height: 100%;
		top:0!important;
		bottom:0!important;
		right:0!important;
		left:0!important;
		margin :0!important;
	   cursor: crosshair !important;
	 background-image: url("./images/grid.png");
	 background-color: #263238;
    }
	


#table-container {
  height: 0%;
  position: relative;
  box-shadow: inset 0 8px 8px -8px #696868;
  display:none;
}



	
	
#overlayers {display:none;}
#toolslayers{display:none;}
#toolsearches{display:none;}
#imageanalysis{display:none;}
#userlayers{display:none;color:black;}
#spatialanalysis{display:none;}



.leaflet-bar a {
      background-color: #f4f4f4;
      display: inline;
    }

    .leaflet-bar a:hover {
      display: inline;
    }
/*elevation profile html*/



.leaflet-map-pane {
    z-index: 2 !important;
}



#labeluserid
{
	
font-size:18px;	

margin-right: 5px;

}




.modal-body {
max-height: 430px!important;
}



.leaflet-routing-alt 

{
	background-color: #505053;
}

</style>


	<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.3/dist/leaflet.css" integrity="sha256-kLaT2GOSpHechhsozzB+flnD+zUyjE2LlfWPgU04xyI=" crossorigin="" />
<script src="https://unpkg.com/leaflet@1.9.3/dist/leaflet.js" integrity="sha256-WBkoXOwTeyKclOHuWtc+i2uENFpDZ9YPdf5Hf+D7ewM=" crossorigin=""></script>

<link rel="stylesheet" href="./css/styles.css"/>
<!--link  rel="stylesheet" href="./css/fontawesome/css/fontawesome.css"/-->
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.2.1/css/all.min.css" integrity="sha512-MV7K8+y+gLIBoVD59lQIYicR65iaqukzvf/nwasF0nqhPay5w/9lJmVM2hMDcnK1OnMGCdVK+iQrJ7lzPJQd1w==" crossorigin="anonymous" referrerpolicy="no-referrer" />

<script  src="./js/turf@6.5.0turf.min.js" crossorigin=""></script>

<script src="./js/aliencraft-map-plugins-leaflet-0.7.js" crossorigin="" ></script>


</head>



	
<body>






<script>




 if (window.location.protocol === "http:") 
 {
location.href = location.href.replace("http://", "https://");
 }

        (function onLoad(){
		
			
          function id(v){ return document.getElementById(v); }
          function loadbar() {
            var ovrl = id("overlay"),
                prog = id("progress"),
                stat = id("progstat"),
                img = document.images,
                c = 0,
                tot = img.length;
            if(tot == 0) return doneLoading();

            function imgLoaded(){
              c += 1;
              var perc = ((100/tot*c) << 0) +"%";
              prog.style.width = perc;
              stat.innerHTML = "Loading "+ perc;
              if(c===tot) return doneLoading();
            }
            function doneLoading(){
              ovrl.style.opacity = 0;
              setTimeout(function(){ 
                ovrl.style.display = "none";
              }, 500);
            }
            for(var i=0; i<tot; i++) {
              var tImg     = new Image();
              tImg.onload  = imgLoaded;
              tImg.onerror = imgLoaded;
              tImg.src     = img[i].src;
            }    
          }
          document.addEventListener('DOMContentLoaded', loadbar, false);
		  

	
        }())
		

    </script>








<div id="sidebars" class="sidebarss" >  
<div class="tab">
			 <a class="btn tablinks" id="identifybtn"   title="identify on map" href="JavaScript:void(0);"  onclick="opentab(event, 'identifybtn')" ><i class="fas fa-info-circle"></i></a>
             <a class="btn tablinks" id="spatialanalysisbtn"   title="spatial analysis" href="JavaScript:void(0);"  onclick="opentab(event, 'spatialanalysis')" ><i class="fa fa-bar-chart"></i></a>
<button class="closetabs" id="closetabs" >??</button>
</div>







<div id="spatialanalysis" class="tabcontent" >
            <hr>
            <h6 class="mdl-color-text--primary mdl-typography--text-left mdl-typography--text-uppercase"><b>Analyze data</b></h6>
            <br>
            <br>
            <br>

			    <div>
               <button onclick="initializenetwork()" class="btn" value="pathnetworkfunction" id="pathnetworkfunction" style="display:none;height: 56px;">Shortest path on Multiline network &emsp;
			   <i class='fas fa-route'></i></button>
            </div>

			   <br>
			    <div>
               <button onclick="choosevoronoilayer()" class="btn" value="voronoifunction" id="voronoifunction" style="display:none;">Create Voronoi Diagram
			   &emsp; <span class="voronoiicon"></span> </button>
            </div>
			
			   <br>
            <div>
               <button onclick="choosepointlayer()" class="btn" value="nearestfunction" id="nearestfunction" style="display:none;"> Closest point from multiPoint
			 &emsp; <span class="shorticon"></span> </i></button>
            </div>
            <br>
        
            <div>
               <button onclick="prepareheatmapdata()" class="btn" value="prepareheatmapdata" id="prepareheatmapdata" style="display:none;">Create Heatmap
              &emsp; <span class="heatmapicon"></span> 	</button>
            </div>
            <br>
      
</div>


<br>



 </div>





<div id="map">

 <div id="overlay">
            <div id="progstat"></div>
            <div id="progress"></div>
       
 </div>       
 
 
 
   <div id="heatmapdivpanel" style="display:none; max-width:300px;">
            Heat Dispersal<input  type="range" id="radiusheat"    min="10"  max="30"  step="0.01" class="mdl-slider mdl-js-slider" />
         </div>
		 
		 
 </div> 






	<div id="identifeatures" style="display:block;">

	<div id="identifyaccordion">
	<h3>Spatial Query <i class="fa fa-cog" aria-hidden="true"></i></h3>

<div class="custom-control custom-checkbox custom-control-inline">



	
  	<div id="identifyswitches">
	
 <br> 
			  
	

  	   	  <label id="identifyoverpassfeatures" class="mdl-card--identifyoverpassfeatures mdl-switch mdl-js-switch mdl-js-ripple-effect" for="identifyoverpassfeaturesid">
                <input type="checkbox" id="identifyoverpassfeaturesid" class="mdl-switch__input" checked >
                <span class="mdl-switch__label">OSM Data</span>
          </label>
			
	 
<select name="myosmvalue" id="overpassquerylayerSelect" style="max-width: 140px;">
	<!--input list="osmvalues" name="myosmvalue" id="overpassquerylayerSelect"/-->
	<!--datalist id="osmvalues"-->
    <option value="none" selected disabled>Select an Option</option>
    <option value="highway" selected >Roads</option>
	<option value="building">building</option>
    <option value="power">Power</option>
	<option value="waterway">waterway</option>
	<option value="tourism=attraction">attraction</option>
	<option value="amenity=parking">parking</option>
	<!--option value="">Input your own value<input id="overpassquerylayerSelectinput" oninput="getvaluefrominput(this.value)"/></option-->
	
<!--/datalist-->
  </select>

<br>

<div class="pad5" id="Rectangle_size" style="display:block">
<input type="range" id="Rectangle_value"  min="2000" max="20000" value="10000" class="slider" oninput="updateTextInput(this.value);" />
<label for="Rectangle_value" id="Rectangle_valuelabel"></label>
</div>
<br>
		   
</div>

</div>

</div><!--//accordion-->


<div id= "countfeatures">Click on map to identify</div>

   <div id="closebuttonidentify">X</div>
	


 <div class="leftpanel" id="leftpanel" style="display:none;">
			
<br>
<div id="toolbar" >
<ul></ul>
</div>
</div><!--end of left panel-->


</div><!--identifeatures-->






















<div id="sidebar">

<div id="mapinfoheader" > 
</div>
<div id="mapinfo">	

</div>


</div>




<style>



	  
	  .btn span.voronoiicon {
    background: url(./images/voronoi-diagram.png) no-repeat;
    float: right;
    width: 40px;
    height: 30px;
}

	  .btn span.heatmapicon {
    background: url(./images/heatmap2.png) no-repeat;
    float: right;
    width: 40px;
    height: 30px;
}


	  .btn span.shorticon {
    background: url(./images/cry.png) no-repeat;
    float: right;
    width: 40px;
    height: 30px;
}

#overlay{
  position:fixed;
  z-index:99999;
  top:0;
  left:0;
  bottom:0;
  right:0;
  background:rgba(0,0,0,0.9);
  transition: 1s 0.4s;
}
#progress{
  height:4px;
  background:#fff;
  position:absolute;
  width:0;
  top:50%;
  transition: 1s;
}
#progstat{
  font-size:2.5em;
  letter-spacing: 3px;
  position:absolute;
  top:50%;
  margin-top:-40px;
  width:100%;
  text-align:center;
  color:#fff;
}


#closetabs{    
width: 39px;
height: 29px;
border-radius: 6px;
border-style: double;
border-color:black;
color: black;
border-width: 1px;
}

#sidebars {

    background-color: rgb(205, 223, 227)!important;
}

#sidebars a, .sidebars a{
    text-decoration: none!important;

}


#sidebars a :hover, .sidebars a :hover
{

		 opacity: 1.0 !important;
        box-shadow:none!important;
						
}





#closebuttonidentify{

z-index: 99;
    position: fixed;
    top: 175px;
    height: 28px!important;
    left: 242px;
	/*right:27px;*/
    width: 33px;
    /* overflow-x: hidden; */
    /* overflow-y: auto; */
    font-size: 22px;
    font-weight: 800;
    text-align: center;
    padding-top: 3px;
    border-radius: 12px;
    background-color: #deac6e;
	   display:none;
}
 





#identifyaccordion

{
position:absolute;
background-color:#d6b6d4;
    z-index: 100;
	top: 105px;
	/*right: 15px;*/
	width: 280px;
	/*background-color: #fff;
	border-radius: 15px;*/
	padding: 10px;
	/*border: solid 4px #3e545f;*/
	font-size: 14px;
   font-weight: 500;
    -moz-box-shadow: 10px;
	border-style: solid;
    border-width: 1px;
	display:none;
	margin-left: 5px;
}








#countfeatures

{
position:fixed;
background-color:#d6b6d4;
    z-index: 90;
	top: 165px;
	height: 45px!important;
	/*right: 15px;*/
	width: 280px;
	/*background-color: #fff;
	border-radius: 15px;*/
	padding: 10px;
	/*border: solid 4px #3e545f;*/
	font-size: 14px;
   font-weight: 500;
    -moz-box-shadow: 10px;
	color: black;
	display:none;
	border-style: double;
	margin-left: 5px;
}





#leftpanel {
z-index: 90;
	position: absolute;
	top: 238px;
	max-height: auto!important;
	/*right: 15px;*/
	width: 280px;
	border-radius: 5px;
	padding: 10px;
	/*border: solid 4px #3e545f;*/
	/*-webkit-box-shadow: 10px 10px 5px 0px rgba(0,0,0,0.5);*/
	/*-moz-box-shadow: 10px 10px 5px 0px rgba(0,0,0,0.5);*/
	/*box-shadow: 10px 10px 5px 0px rgba(0,0,0,0.5);*/
	max-height: 30%;
	overflow-x: hidden;
	overflow-y: auto;
	font-size: 14px;
   font-weight: 500;
    -moz-box-shadow: 10px;
	font-family: -apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica Neue,Arial,Noto Sans,sans-serif;
    /*background-color: #f0fdffde;*/
	background-color: #d6b6d4;
    color: black;
	margin-left: 5px;
}


.leaflet-routing-container h3 {
    font-size: 25px!important;
    font-weight: normal!important;
}


.leaflet-control-layers
{
	
	background-color: #40f121!important;
	width: 200px;
    font-size: 16px;
}



</style>


</body>


<script  src="./js/index.js"></script>


</html>
