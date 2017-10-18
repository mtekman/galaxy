"use strict";function _interopRequireDefault(e){return e&&e.__esModule?e:{default:e}}Object.defineProperty(exports,"__esModule",{value:!0});var _localization=require("utils/localization"),_localization2=_interopRequireDefault(_localization),_historyModel=require("mvc/history/history-model"),_historyModel2=_interopRequireDefault(_historyModel),_historyViewEdit=require("mvc/history/history-view-edit"),_historyViewEdit2=_interopRequireDefault(_historyViewEdit),_copyDialog=require("mvc/history/copy-dialog"),_copyDialog2=_interopRequireDefault(_copyDialog),_errorModal=require("mvc/ui/error-modal"),_errorModal2=_interopRequireDefault(_errorModal),_baseMvc=require("mvc/base-mvc"),_baseMvc2=_interopRequireDefault(_baseMvc),_ajaxQueue=require("utils/ajax-queue"),_ajaxQueue2=_interopRequireDefault(_ajaxQueue);require("ui/mode-button"),require("ui/search-input");var logNamespace="history",HistoryViewColumn=Backbone.View.extend(_baseMvc2.default.LoggableMixin).extend({_logNamespace:logNamespace,tagName:"div",className:"history-column flex-column flex-row-container",id:function(){return this.model?"history-column-"+this.model.get("id"):""},initialize:function(e){e=e||{},this.purgeAllowed=!_.isUndefined(e.purgeAllowed)&&e.purgeAllowed,this.panel=e.panel||this.createPanel(e),this.setUpListeners()},createPanel:function(e){return new _historyViewEdit2.default.HistoryViewEdit(_.defaults(e,{model:this.model,purgeAllowed:this.purgeAllowed,dragItems:!0,$scrollContainer:function(){return this.$el}}))},setUpListeners:function(){var e=this;this.once("rendered",function(){e.trigger("rendered:initial",e)}),this.setUpPanelListeners()},setUpPanelListeners:function(){var e=this;this.listenTo(this.panel,{rendered:function(){e.trigger("rendered",e)},"view:expanded view:rendered":function(e){e.$(".rerun-btn").off()}},this)},inView:function(e,t){var i=this.$el.offset().left;return!(i+this.$el.width()<e)&&!(i>t)},$panel:function(){return this.$(".history-panel")},render:function(e){e=void 0!==e?e:"fast";var t=this.model?this.model.toJSON():{};return this.$el.html(this.template(t)),this.renderPanel(e),this.panel.$el.css("display","flex"),this.setUpBehaviors(),this},setUpBehaviors:function(){},template:function(e){return e=_.extend(e||{},{isCurrentHistory:this.currentHistory}),$(['<div class="panel-controls clear flex-row">',this.controlsLeftTemplate({history:e,view:this}),this.controlsRightTemplate({history:e,view:this}),"</div>",'<div class="inner flex-row flex-column-container">','<div id="history-',e.id,'" class="history-column history-panel flex-column"></div>',"</div>"].join(""))},renderPanel:function(e){return e=void 0!==e?e:"fast",this.panel.setElement(this.$panel()).render(e),this.currentHistory&&this.panel.$list().before(this.panel._renderDropTargetHelp()),this},events:{"click .switch-to.btn":function(){this.model.setAsCurrent()},"click .delete-history":function(){var e=this;this.model._delete().done(function(t){e.render()})},"click .undelete-history":function(){var e=this;this.model.undelete().done(function(t){e.render()})},"click .purge-history":function(){if(confirm((0,_localization2.default)("This will permanently remove the data. Are you sure?"))){var e=this;this.model.purge().done(function(t){e.render()})}},"click .copy-history":"copy"},copy:function(){(0,_copyDialog2.default)(this.model)},controlsLeftTemplate:_.template(['<div class="pull-left">',"<% if( data.history.isCurrentHistory ){ %>",'<strong class="current-label">',(0,_localization2.default)("Current History"),"</strong>","<% } else { %>",'<button class="switch-to btn btn-default">',(0,_localization2.default)("Switch to"),"</button>","<% } %>","</div>"].join(""),{variable:"data"}),controlsRightTemplate:_.template(['<div class="pull-right">',"<% if( !data.history.purged ){ %>",'<div class="panel-menu btn-group">','<button type="button" class="btn btn-default dropdown-toggle" data-toggle="dropdown">','<span class="caret"></span>',"</button>",'<ul class="dropdown-menu pull-right" role="menu">',"<% if( !data.history.deleted ){ %>",'<li><a href="javascript:void(0);" class="copy-history">',(0,_localization2.default)("Copy"),"</a></li>",'<li><a href="javascript:void(0);" class="delete-history">',(0,_localization2.default)("Delete"),"</a></li>","<% } else /* if is deleted */ { %>",'<li><a href="javascript:void(0);" class="undelete-history">',(0,_localization2.default)("Undelete"),"</a></li>","<% } %>","<% if( data.view.purgeAllowed ){ %>",'<li><a href="javascript:void(0);" class="purge-history">',(0,_localization2.default)("Purge"),"</a></li>","<% } %>","</ul>","</div>","<% } %>","</div>"].join(""),{variable:"data"}),toString:function(){return"HistoryViewColumn("+(this.panel?this.panel:"")+")"}}),MultiPanelColumns=Backbone.View.extend(_baseMvc2.default.LoggableMixin).extend({_logNamespace:logNamespace,className:"multi-panel-history",initialize:function(e){e=e||{},this.log(this+".init",e),this.$el.addClass(this.className),this.options={columnWidth:312,borderWidth:1,columnGap:8,headerHeight:29,footerHeight:0,controlsHeight:20},this.perPage=e.perPage||10,this.hdaQueue=new _ajaxQueue2.default.NamedAjaxQueue([],!1),this.collection=null,this.columnMap={},this.columnOptions=e.columnOptions||{},this.historySearch=null,this.datasetSearch=null,this.setCollection(e.histories),this.setUpListeners()},setUpListeners:function(){var e=this;this.on("end-of-scroll",function(){e.collection.fetchMore()})},setCollection:function(e){return this.stopListening(this.collection),this.collection=e||new _historyModel2.default.HistoryCollection,this.setUpCollectionListeners(),this.createColumns(),this.hdaQueue.clear(),this.trigger("new-collection",this),this},addModels:function(e,t,i){i=i||{};var o=this;return(e=_.isArray(e)?e:[e]).forEach(function(e){o.addColumn(e,!1)}),this},setUpCollectionListeners:function(){var e=this;e.listenTo(e.collection,{error:e.errorHandler,add:e.addModels,"all-fetched":e._postFetchAll,"new-current":e.addAsCurrentColumn,"set-as-current":e.setCurrentHistory,"change:deleted change:purged":e.handleDeletedHistory,sort:function(){e.renderColumns(0)}})},_postFetchAll:function(e){if(this.$(".histories-loading-indicator").remove(),!this.historySearch){var t=this.$(".outer-middle");t.scrollLeft(t.scrollLeft()+24)}},setCurrentHistory:function(e){this.log("setCurrentHistory:",e);var t=_.findWhere(this.columnMap,{currentHistory:!0});t&&(t.currentHistory=!1,t.$el.height(""));var i=this.columnMap[this.collection.currentHistoryId];return i.currentHistory=!0,this.collection.sort(),this._recalcFirstColumnHeight(),i},handleDeletedHistory:function(e){if(e.get("deleted")||e.get("purged")){this.log("handleDeletedHistory",this.collection.includeDeleted,e);var t=this,i=t.columnMap[e.id];if(!i)return;i.model.id===this.collection.currentHistoryId||t.collection.includeDeleted||t.removeColumn(i)}},errorHandler:function(e,t,i){if(!t||0!==t.status||0!==t.readyState){if(this.error(e,t,i),_.isString(e)&&_.isString(t)){var o=e,n=t;return _errorModal2.default.errorModal(o,n,i)}return t&&502===t.status?_errorModal2.default.badGatewayErrorModal():_errorModal2.default.ajaxErrorModal(e,t,i)}},_ajaxErrorHandler:function(){_errorModal2.default.ajaxErrorModal.apply(null,_.toArray(arguments))},create:function(e){return this.collection.create({current:!0})},createColumns:function(e,t){t=t||this.options.columnOptions;var i=this;i.columnMap={},i.collection.each(function(e,o){var n=i.createColumn(e,t);i.columnMap[e.id]=n})},createColumn:function(e,t){t=_.extend({},t,{model:e,purgeAllowed:Galaxy.config.allow_user_dataset_purge});var i=new HistoryViewColumn(t);return e.id===this.collection.currentHistoryId&&(i.currentHistory=!0),this.setUpColumnListeners(i),this.datasetSearch&&(i.panel.searchItems(this.datasetSearch),this.queueHdaFetchDetails(i)),i},addColumn:function(e,t){t=void 0===t||t;var i=this.createColumn(e);return this.columnMap[e.id]=i,t&&this.renderColumns(),i},addAsCurrentColumn:function(e,t,i){var o=this,n=this.addColumn(e,!1);return this.setCurrentHistory(e),n.once("rendered",function(){o.queueHdaFetch(n)}),n},removeColumn:function(e,t){if(t=void 0===t||t,this.log("removeColumn",e),e){var i=this,o=this.options.columnWidth+this.options.columnGap;e.$el.fadeOut("fast",function(){t&&($(this).remove(),i.$(".middle").width(i.$(".middle").width()-o),i.checkColumnsInView(),i._recalcFirstColumnHeight()),i.stopListening(e.panel),i.stopListening(e),delete i.columnMap[e.model.id],e.remove()})}},setUpColumnListeners:function(e){var t=this;t.listenTo(e,{"in-view":t.queueHdaFetch}),t.listenTo(e.panel,{"view:draggable:dragstart":function(e,i,o,n){t._dropData=JSON.parse(e.dataTransfer.getData("text")),t.currentColumnDropTargetOn()},"view:draggable:dragend":function(e,i,o,n){t._dropData=null,t.currentColumnDropTargetOff()},"droptarget:drop":function(e,i,o){var n=t._dropData.filter(function(e){return o.model.contents.isCopyable(e)});t._dropData=null;var r=new _ajaxQueue2.default.NamedAjaxQueue;0!==o.model.contents.currentPage&&r.add({name:"fetch-front-page",fn:function(){return o.model.contents.fetchPage(0)}}),n.reverse().forEach(function(e){r.add({name:"copy-"+e.id,fn:function(){return o.model.contents.copy(e)}})}),r.start(),r.done(function(e){o.model.fetch()})}})},columnMapLength:function(){return Object.keys(this.columnMap).length},sortedFilteredColumns:function(e){return(e=e||this.filters)&&e.length?this.sortedColumns().filter(function(t,i){return t.currentHistory||_.every(e.map(function(e){return e.call(t)}))}):this.sortedColumns()},sortedColumns:function(){var e=this;return this.collection.map(function(t,i){return e.columnMap[t.id]})},render:function(e){e=void 0!==e?e:this.fxSpeed;var t=this;return t.log(t+".render"),t.$el.html(t.mainTemplate(t)),t.renderColumns(e),t.setUpBehaviors(),t.trigger("rendered",t),t},renderColumns:function(e){e=_.isNumber(e)?e:this.fxSpeed;var t=this,i=t.sortedFilteredColumns(),o=t.$(".middle").empty();return t._addColumns(i,e),t.collection.allFetched||o.append(t.loadingIndicatorTemplate(t)),t.trigger("columns-rendered",i,t),t.datasetSearch&&i.length<=1||(t.checkColumnsInView(),t._recalcFirstColumnHeight()),i},_addColumns:function(e,t){t=_.isNumber(t)?t:this.fxSpeed;var i=this.$(".middle"),o=i.children(".history-column").length;i.width(this._calcMiddleWidth(e.length+o)),e.forEach(function(e,o){e.delegateEvents().render(t).$el.appendTo(i)})},_calcMiddleWidth:function(e){return e*(this.options.columnWidth+this.options.columnGap)+this.options.columnGap+16},queueHdaFetch:function(e){var t=e.model.contents;if(0===t.length&&e.model.contentsShown()){var i={silent:!0},o=_.values(t.storage.allExpanded()).join();o&&(i.details=o),this.hdaQueue.add({name:e.model.id,fn:function(){return t.fetchCurrentPage(i).done(function(){e.panel.renderItems()})}}),this.hdaQueue.running||this.hdaQueue.start()}},queueHdaFetchDetails:function(e){var t=e.model.contents;!(0===t.length&&e.model.contentsShown())&&t.haveDetails()||(this.hdaQueue.add({name:e.model.id,fn:function(){return t.progressivelyFetchDetails().done(function(){e.panel._renderEmptyMessage()})}}),this.hdaQueue.running||this.hdaQueue.start())},renderInfo:function(e){return this.$(".header .header-info").text(e)},events:{"click .done.btn":"close","click .create-new.btn":"create","click #include-deleted":"_clickToggleDeletedHistories","click .order .set-order":"_chooseOrder","click #toggle-deleted":"_clickToggleDeletedDatasets","click #toggle-hidden":"_clickToggleHiddenDatasets"},close:function(e){window.location=Galaxy.root},_clickToggleDeletedHistories:function(e){this.toggleDeletedHistories($(e.currentTarget).is(":checked")),this.toggleOptionsPopover()},toggleDeletedHistories:function(e){window.location=e?Galaxy.root+"history/view_multiple?include_deleted_histories=True":Galaxy.root+"history/view_multiple"},_clickToggleDeletedDatasets:function(e){this.toggleDeletedDatasets($(e.currentTarget).is(":checked")),this.toggleOptionsPopover()},toggleDeletedDatasets:function(e){e=void 0!==e&&e,this.sortedFilteredColumns().forEach(function(t,i){_.delay(function(){t.panel.toggleShowDeleted(e,!1)},200*i)})},_clickToggleHiddenDatasets:function(e){this.toggleHiddenDatasets($(e.currentTarget).is(":checked")),this.toggleOptionsPopover()},toggleHiddenDatasets:function(e){e=void 0!==e&&e,this.sortedFilteredColumns().forEach(function(t,i){_.delay(function(){t.panel.toggleShowHidden(e,!1)},200*i)})},_chooseOrder:function(e){var t=this,i=t.collection,o=$(e.currentTarget).data("order");t.$(".current-order").text(t.orderDescriptions[o]),t.toggleOptionsPopover(),i.setOrder(o);var n=i.slice(0,1);i.fetchFirst().done(function(){i.unshift(n,{silent:!0}),t.createColumns(),t.hdaQueue.clear(),t.render()}),t.once("columns-rendered",t._scrollLeft)},_scrollLeft:function(e){e=_.isNumber(e)?e:0,this.$(".outer-middle").scrollLeft(e)},setUpBehaviors:function(){var e=this;e._moreOptionsPopover(),e.$("#search-histories").searchInput({name:"search-histories",placeholder:(0,_localization2.default)("search histories"),onfirstsearch:function(t){e.$("#search-histories").searchInput("toggle-loading"),e.renderInfo((0,_localization2.default)("loading all histories for search")),e.collection.fetchAll().done(function(){e.$("#search-histories").searchInput("toggle-loading"),e.renderInfo("")})},onsearch:function(t){e.historySearch=t,e.filters=[function(){return this.model.matchesAll(e.historySearch)}],e.renderColumns(0)},onclear:function(t){e.historySearch=null,e.filters=[],e.renderColumns(0)}}),e.$("#search-datasets").searchInput({name:"search-datasets",placeholder:(0,_localization2.default)("search all datasets"),onfirstsearch:function(t){e.hdaQueue.clear(),e.$("#search-datasets").searchInput("toggle-loading"),e.datasetSearch=t,e.sortedFilteredColumns().forEach(function(i){i.panel.searchItems(t),e.queueHdaFetchDetails(i)}),e.hdaQueue.progress(function(t){e.renderInfo([(0,_localization2.default)("searching"),t.curr+1,(0,_localization2.default)("of"),t.total].join(" "))}),e.hdaQueue.deferred.done(function(){e.renderInfo(""),e.$("#search-datasets").searchInput("toggle-loading")})},onsearch:function(t){e.datasetSearch=t,e.sortedFilteredColumns().forEach(function(e){e.panel.searchItems(t)})},onclear:function(t){e.datasetSearch=null,e.sortedFilteredColumns().forEach(function(e){e.panel.clearSearch()})}}),$(window).resize(function(){e._recalcFirstColumnHeight()});var t=_.debounce(function(){var t=e._viewport();e.checkColumnsInView(t),e.checkForEndOfScroll(t)},100);this.$(".middle").parent().scroll(t)},_moreOptionsPopover:function(){return this.$(".open-more-options.btn").popover({container:".header",placement:"bottom",html:!0,content:$(this.optionsPopoverTemplate(this))})},toggleOptionsPopover:function(e){this.$(".open-more-options.btn").popover("toggle")},_recalcFirstColumnHeight:function(){var e=this.$(".history-column").first(),t=this.$(".middle").height(),i=e.find(".panel-controls").height();e.height(t).find(".inner").height(t-i)},_viewport:function(){var e=this.$(".middle").parent(),t=e.offset().left;return{left:t,right:t+e.width()}},columnsInView:function(e){var t=e||this._viewport();return this.sortedFilteredColumns().filter(function(e){return e.currentHistory||e.inView(t.left,t.right)})},checkColumnsInView:function(){this.columnsInView().forEach(function(e){e.trigger("in-view",e)})},checkForEndOfScroll:function(e){e=e||this._viewport();var t=this.$(".middle");t.parent().scrollLeft()+e.right>=t.width()-16&&this.trigger("end-of-scroll")},currentColumnDropTargetOn:function(){var e=this.columnMap[this.collection.currentHistoryId];e&&(e.panel.dataDropped=function(e){},e.panel.dropTargetOn())},currentColumnDropTargetOff:function(){var e=this.columnMap[this.collection.currentHistoryId];e&&(e.panel.dataDropped=_historyViewEdit2.default.HistoryViewEdit.prototype.dataDrop,e.panel.dropTarget=!1,e.panel.$(".history-drop-target").remove())},toString:function(){return"MultiPanelColumns("+(this.columns?this.columns.length:0)+")"},mainTemplate:_.template(['<div class="header flex-column-container">','<div class="control-column control-column-left flex-column">','<button class="done btn btn-default" tabindex="1">',(0,_localization2.default)("Done"),"</button>",'<div id="search-histories" class="search-control"></div>','<div id="search-datasets" class="search-control"></div>','<a class="open-more-options btn btn-default" tabindex="3">','<span class="fa fa-ellipsis-h"></span>',"</a>","</div>",'<div class="control-column control-column-center flex-column">','<div class="header-info">',"</div>","</div>",'<div class="control-column control-column-right flex-column">','<button class="create-new btn btn-default" tabindex="4">',(0,_localization2.default)("Create new"),"</button> ","</div>","</div>",'<div class="outer-middle flex-row flex-row-container">','<div class="middle flex-column-container flex-row"></div>',"</div>",'<div class="footer flex-column-container"></div>'].join(""),{variable:"view"}),loadingIndicatorTemplate:_.template(['<div class="histories-loading-indicator">','<span class="fa fa-spin fa-spinner"></span>',(0,_localization2.default)("Loading histories"),"...","</div>"].join(""),{variable:"view"}),orderDescriptions:{update_time:(0,_localization2.default)("most recent first"),"update_time-asc":(0,_localization2.default)("least recent first"),name:(0,_localization2.default)("name, a to z"),"name-dsc":(0,_localization2.default)("name, z to a"),size:(0,_localization2.default)("size, large to small"),"size-asc":(0,_localization2.default)("size, small to large")},optionsPopoverTemplate:_.template(['<div class="more-options">','<div class="order btn-group">','<button type="button" class="btn btn-default dropdown-toggle" data-toggle="dropdown">',(0,_localization2.default)("Order histories by")+" ",'<span class="current-order"><%- view.orderDescriptions[ view.collection.order ] %></span> ','<span class="caret"></span>',"</button>",'<ul class="dropdown-menu" role="menu">',"<% _.each( view.orderDescriptions, function( text, order ){ %>",'<li><a href="javascript:void(0);" class="set-order" data-order="<%- order %>">',"<%- text %>","</a></li>","<% }); %>","</ul>","</div>",'<div class="checkbox"><label><input id="include-deleted" type="checkbox"','<%= view.collection.includeDeleted? " checked" : "" %>>',(0,_localization2.default)("Include deleted histories"),"</label></div>","<hr />",'<div class="checkbox"><label><input id="toggle-deleted" type="checkbox">',(0,_localization2.default)("Include deleted datasets"),"</label></div>",'<div class="checkbox"><label><input id="toggle-hidden" type="checkbox">',(0,_localization2.default)("Include hidden datasets"),"</label></div>","</div>"].join(""),{variable:"view"})});exports.default={MultiPanelColumns:MultiPanelColumns};
//# sourceMappingURL=../../../maps/mvc/history/multi-panel.js.map
