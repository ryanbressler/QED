var View = require('./view');
var template = require('./templates/branch');
var BranchMatrix = require('../models/branchMatrix');

module.exports = View.extend({

  model:BranchMatrix,
  template:template,

  initialize : function() {
  	_.bindAll(this, 'getRenderData', 'afterRender', 'renderGrid');
  },
  
  getRenderData : function() {},

  afterRender: function() {
  	var _this = this;
    this.$el.addClass('row-fluid');
    //this.model.bind('load',_this.renderGrid);
    this.model.on('load',_this.renderGrid);
  },

  featureScat : function(feature){
    console.log(feature);
  },

  renderGrid : function(){
    var me = this;
    var scatdata = this.model.getD3Data();
    var length = scatdata.length;
    var ignore_keys = ['label','type','source','feature_id','nByi',"feature","featureID"];
    var height = 200;
    var width = 1000;
    var x = function(d) {
          return d.r1*(width/length);
      };
    var y = function(d) {
          return d.r2*(height/length);
      };
   
    var colors = ["#D7191C", "#FDAE61", "#CCCCAC", "#ABD9E9", "#2C7BB6"];
    var color = function(d){
        if ( ! colors.hasOwnProperty(parseInt(d.termcat))){
          console.log("bad color "+ (d.termcat-1));
          return "#000000";
        }
        return colors[d.termcat-1];
      };

    var clearstate = function ()
    {
      $(".case-name").html("");
        $(".branch-output").html("");
        $(".pathway-output").html("");
        $(".admix-output").html("");
        me.pc.color("#447");
        me.pc.clear("highlight");
        me.pc.render();
    }
    var showpathways=function(data){
      if (data.length==0){
        $(".pathway-output").html("No Pathways Found.");
        return;
      }
      $(".pathway-output").html("<div style='font-size: 14px;'>Pathways Containing "+data[0].GENE+"</div>");
      var grid = d3.select(".pathway-output");
      grid.selectAll("div")
        .data(data)
        .enter()
        .append("div")
        .text(function(d){ return d.PATHWAY;});

    };

    var highlightf = function(d){me.pc.highlight(me.model.filterFeatures([d[0]]));};

    var admixcolor = d3.scale.linear()
      .domain([0, 1])
      .range(["steelblue", "brown"])
      .interpolate(d3.interpolateLab);

    var loadAdMix = function(name){
      d3.tsv("/svc/data/analysis/admix4.txt?cols="+name,function(data){
        var grid = d3.select(".admix-output");
        grid.selectAll("span")
          .data(data)
          .enter()
          .append("span")
          .text(function(d){ return d.SAMPLE+":"+d[name]+" ";})
          .style("color",function(d){ return admixcolor(d[name]); });

      });
    };

    var showFvsT = function (selector,data) {
      $(selector).html("");
      var grid = d3.select(selector);
        grid.selectAll("div")
          .data(data)
          .enter()
          .append("div")
          .text(function(d){ return d[1];})
          .style("width",function(d){return d[1]+"px";})
          .style("float","left")
          .style("background-color",function(d){ console.log(d[0]); return color({termcat:d[0]}); });
    };

    var loadFvsTermCat = function(name,target,selector){
      var fmsvcbase = "";
      var dataset_id = me.model.get("dataset_id");
      var bestlength = 0;
      d3.json("/svc/data/domains/feature_matrices",function(fms){
        for (var i = fms.files.length - 1; i >= 0; i--) {
          if (dataset_id.indexOf(fms.files[i].label) === 0 && fms.files[i].label.length > bestlength ) {
            fmsvcbase=fms.files[i].uri;
            bestlength=fms.files[i].label.length;
          }
        }

        d3.tsv("/svc"+fmsvcbase+"?rows=B:MRGE:Strict_Hypertension_Related:NB::::,"+target+","+name,function(data){
          var truevs = [];
          var falsevs = [];
          var pcdata = me.model.get('branches');
          var countByCase={};
          var maxCount = 0;
          var scaterdata = [];
          iByn={};

          for (var i = 0; i < data.length; i++) {
              iByn[data[i]["."]]=i;
            }

          var termi=iByn[target];
          var hypei=iByn["B:MRGE:Strict_Hypertension_Related:NB::::"];
          var namei=iByn[name];
          for (var i = pcdata.length - 1; i >= 0; i--) {
            if (pcdata[i][0]==name){
              for (var j = 1; j < pcdata[i].length; j++) {
                var caseid = pcdata[0][j];
                
                var parseBin = function(b) {
                  var v = b.toLowerCase()=="true" ? 1 : 0;
                  return v +0.5*Math.random()-0.25;
                }

                var x = parseFloat(data[namei][caseid]);
                if (name[0]=="B") {
                  x = parseBin(data[namei][caseid]);
                }
                var y = parseFloat(data[termi][caseid]);
                if (target[0]=="B") {
                  y = parseBin(data[termi][caseid]);
                }
                
                scaterdata.push({x:x,
                  y:y,
                  hype:data[hypei][caseid].toLowerCase()=="true",
                  count:pcdata[i][j],
                  caseid:caseid});
                countByCase[pcdata[0][j]]=pcdata[i][j];
                if (pcdata[i][j]>maxCount) {
                  //should this be the global max instead of just for this feature?
                  maxCount = pcdata[i][j];
                }
              }
            }
          }

          
          
          var width = 1000;
          var height = 600; 
          
          var keys = _.difference(Object.keys(data[0]),ignore_keys);
          $(selector).html("");
          var svg = d3.select(selector)
            .append("svg")
            .attr("width", width)
            .attr("height", height);

          var decorate = $("#decorate_fplot").is(':checked')
          
          var opac = function(d) {
            return 0.2+0.4*d.count/maxCount;
          };

          var size = function(d) {
            if (decorate){
              return 4+6*d.count/maxCount;
            }
            return 6;
          };

          var xScale = d3.scale.linear()
                     .domain([d3.min(scaterdata, function(d){return d.x;}), d3.max(scaterdata, function(d){return d.x;})])
                     .range([80, width-40]);
          var yScale = d3.scale.linear()
                     .domain([d3.min(scaterdata, function(d){return d.y;})-0.3, d3.max(scaterdata, function(d){return d.y;})+0.3])
                     .range([40, height-80]);
          var xAxis = d3.svg.axis().orient("bottom").scale(xScale).ticks(4);
          var yAxis = d3.svg.axis().orient("right").scale(yScale).ticks(4);

          var color = function (d) {
            if (decorate) {
              return d.hype ? "red" : "blue";
            }
            return "green";
          }

          // Add the x-axis.
          var xaxiistrans = height -40; 
          svg.append("g")
              .attr("class", "x axis")
              .attr("transform", "translate(0," + xaxiistrans + ")")
              .call(xAxis);

          // Add the y-axis.
          svg.append("g")
              .attr("class", "y axis")
              .attr("transform", "translate(20,0)")
              .call(yAxis);

          // Add an x-axis label.
          svg.append("text")
              .attr("class", "x label")
              .attr("text-anchor", "end")
              .attr("x", width)
              .attr("y", height - 50)
              .text(name);

          // Add a y-axis label.
          svg.append("text")
              .attr("class", "y label")
              .attr("text-anchor", "end")
              .attr("y", 60)
              .attr("x", -20)
              .attr("dy", ".75em")
              .attr("transform", "rotate(-90)")
              .text(target);

          svg.selectAll("circle")
            .data(scaterdata)
            .enter()
            .append("circle")
            .attr("cx", function(d){return xScale(d.x);})
            .attr("cy", function(d){
              if (target != "N:CLIN:Gestational_Age_at_Delivery:NB::::"){
                return yScale(d.y+0.5*Math.random()-0.25);
              }
              return yScale(d.y);
              })
            //.attr("r",10)
            .attr("r", size)
            .attr("fill",color)
            .attr("stroke",color)
            .style('stroke-opacity', .5)
            .style('fill-opacity', .5)
            .append("svg:title")
            .text(function(d) { return d.caseid+" termcat: "+d.y; });
          

        });
      });
    };

    var svg = d3.select(".scat-container")
      .append("svg")
      .attr("width", width)
      .attr("height", height)
      .on("click", clearstate);
      //.margin({ top: 120, left: 80, bottom: 80, right: 80 });
      // Add an x-axis label.
      
    
    svg.selectAll("circle")
      .data(scatdata)
      .enter()
      .append("circle")
      .attr("cx", x)
      .attr("cy", y)
      .attr("r", 10)
      .attr("fill",color)
      .attr("stroke",color)
      .style('stroke-opacity', 0.8)
      .style('fill-opacity', 0.8)
      .on("click", function(d,i){
          clearstate();
          loadAdMix(d.n);
          var features = me.model.getTopFeaturs(d.n);
          var blue_to_brown = d3.scale.pow()
            .exponent(.1)
            .domain([features[features.length-1][1], features[0][1]])
            .range(["#002", "red"])
            .interpolate(d3.interpolateLab);
          me.pc.color(function(d){ return blue_to_brown(d[i]); });
          $(".case-name").html("Admixture and Top Features for "+ d.n);
          var output = d.n + " : ";
          var highlight = [];
          for (var j = 0; j < features.length; j++) {
            highlight.push([features[j][0],features[j][1]]);
            
          }
          //me.pc.highlight(me.model.filterFeatures(highlight));
          var grid = d3.select(".branch-output");

          
          
          grid.selectAll("div")
            .data(highlight)
            .enter()
            .append("div")
            .text(function(d){return d[0]+", "+d[1];})
            .style("color",function(d){ return blue_to_brown(d[1]); })
            .on("click", function(d){
              $(".pathway-output").html("");
              highlightf(d);
              
              loadFvsTermCat(d[0],"N:CLIN:TermCategory:NB::::",".feature-container");
              loadFvsTermCat(d[0],"N:CLIN:Gestational_Age_at_Delivery:NB::::",".feature-container-cont");
              loadFvsTermCat(d[0],"B:CLIN:Preterm:NB::::",".feature-container-binary");


              d3.tsv("/svc/data/analysis/genesets/genesets?rows="+d[0].split(":")[2], showpathways);
              })
            .on("mouseover",highlightf)
            .on("mouseout",function(){me.pc.clear("highlight");});

          //$(".branch-output").html(output);
          me.pc.render();
          d3.event.stopPropagation();
        })
      .append("svg:title")
      .text(function(d) { return d.n; });
    /*svg.append("text")
          .attr("class", "x label")
          .attr("text-anchor", "end")
          .attr("x", width)
          .attr("y", height - 10)
          .text("Case Index in 1st Eigenvector");*/

    // Add a y-axis label.
    svg.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "end")
        .attr("y", 20)
        .attr("x", -20)
        .attr("dy", ".75em")
        .attr("transform", "rotate(-90)")
        .text("Index in 2nd Eigenvector");

    /*svg.selectAll("text")
      .data(scatdata)
      .enter()
      .append("text")
      .text(function(d) {
        return d.n;
      })
      .attr("x", x)
      .attr("y", y)
      .attr("font-family", "sans-serif")
      .attr("font-size", "10px")
      .style('stroke-opacity',.8)
      .style('fill-opacity',.8)
      .attr("fill", "grey");*/


    var pcdata = this.model.get('branches');
    var fnames = [];

    for (var i = 0; i < pcdata.length; i++) {
      fnames.push(pcdata[i][0]);

    }


    $(".feature-search").typeahead({
      source:fnames,
      updater:function (item) {
          highlightf([item]);
          showpathways([items]);
          return item;
      }});

    
    
    
    var keys = _.difference(Object.keys(pcdata[0]),ignore_keys);

    me.pc = d3.parcoords()(".pc-container");
    
    me.pc.dimensions(keys)
      me.pc.data(pcdata)

      //.render()
      .color("#447")
      //.alpha(0.8)
      .margin({ top: 0, left: 0, bottom: 0, right: 0 })
    //me.pc.xscale=d3.scale.linear().domain([0,pcdata[0].length]).range([0,1000]);
    me.pc.render()/*
      .reorderable()
      .brushable()
      .on('brush', function(data){
      me.model.filterNodes(data);
      })*/;
var pcsvg = d3.select(".pc-container")
      .append("svg")
      .attr("width", width)
      .attr("height", 200)
      //.margin({ top: 120, left: 80, bottom: 80, right: 80 });
      // Add an x-axis label.
    pcsvg.append("text")
          .attr("class", "x label")
          .attr("text-anchor", "middle")
          .attr("x", width/2)
          .attr("y", 12)
          .text("Case Index in 1st Eigenvector");
    // Add a y-axis label.
    pcsvg.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "end")
        .attr("y", 20)
        .attr("x", -40)
        .attr("dy", ".75em")
        .attr("transform", "rotate(-90)")
        .text("Local Feature Importance");
		
	}

});