var View = require('./view');
var template = require('./templates/oncovis');
require("../vis/oncovis");
require("../vis/oncovisSlider");

module.exports = View.extend({
    template:template,

    getRenderData:function () {
    },

    afterRender:function () {
        this.$el.addClass('row');
    },

    renderGraph:function () {
        console.log("renderGraph:start");

        var color_fn = function(val) {
            return "blank";
        };

        var optns = {
            bar_width:4,
            column_spacing:1,
            plot_width:3000,
            label_width:70,
            highlight_cls:"mutation",
            color_fn:color_fn,
            columns_by_cluster:oncovisData.columns_by_cluster,
            cluster_labels:oncovisData.cluster_labels,
            row_labels:oncovisData.row_labels
        };

        _.extend(optns, { "bar_height": this.$el.find(".slider_barheight").oncovis_range_value() });
        _.extend(optns, { "row_spacing": this.$el.find(".slider_rowspacing").oncovis_range_value() });
        _.extend(optns, { "bar_width": this.$el.find(".slider_barwidth").oncovis_range_value() });
        _.extend(optns, { "column_spacing": this.$el.find(".slider_barspacing").oncovis_range_value() });
        _.extend(optns, { "cluster_spacing": this.$el.find(".slider_clusterspacing").oncovis_range_value() });
        _.extend(optns, { "label_fontsize": this.$el.find(".slider_fontsize").oncovis_range_value() });

        this.$el.find(".oncovis-container").oncovis(oncovisData.data, optns);

        console.log("renderGraph:end");
    },

    initControls: function() {
        console.log("initControls:start");
        var me = this;
        var visrangeFn = function(property) {
            return function(value) {
                var dim = {};
                dim[property] = value;
                me.$el.find(".oncovis-container").update(dim);
            }
        };

        this.$el.find(".slider_barheight").oncovis_range({ storageId: "slider_barheight", min: 10, max: 50, initialStep: 20, slide: visrangeFn("bar_height") });
        this.$el.find(".slider_rowspacing").oncovis_range({ storageId: "slider_rowspacing", min: 0, max: 50, initialStep: 10, slide: visrangeFn("row_spacing") });
        this.$el.find(".slider_barwidth").oncovis_range({ storageId: "slider_barwidth", min: 1, max: 10, initialStep: 5, slide: visrangeFn("bar_width") });
        this.$el.find(".slider_barspacing").oncovis_range({ storageId: "slider_barspacing", min: 0, max: 10, initialStep: 2, slide: visrangeFn("column_spacing") });
        this.$el.find(".slider_clusterspacing").oncovis_range({ storageId: "slider_clusterspacing", min: 0, max: 50, initialStep: 10, slide: visrangeFn("cluster_spacing") });
        this.$el.find(".slider_fontsize").oncovis_range({ storageId: "slider_fontsize", min: 5, max: 21, initialStep: 14, slide: visrangeFn("label_fontsize") });

        this.$el.find(".reset-sliders").click(function() {
            me.$el.find(".slider_barheight").oncovis_range_reset();
            me.$el.find(".slider_rowspacing").oncovis_range_reset();
            me.$el.find(".slider_barwidth").oncovis_range_reset();
            me.$el.find(".slider_barspacing").oncovis_range_reset();
            me.$el.find(".slider_clusterspacing").oncovis_range_reset();
            me.$el.find(".slider_fontsize").oncovis_range_reset();
        });

        console.log("initControls:end");
    }

});