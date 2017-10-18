"use strict";function _interopRequireDefault(e){return e&&e.__esModule?e:{default:e}}Object.defineProperty(exports,"__esModule",{value:!0});var _historyView=require("mvc/history/history-view"),_historyView2=_interopRequireDefault(_historyView),_hdaLi=require("mvc/history/hda-li"),_hdaLi2=_interopRequireDefault(_hdaLi),_hdcaLi=require("mvc/history/hdca-li"),_hdcaLi2=_interopRequireDefault(_hdcaLi),_baseMvc=require("mvc/base-mvc"),_baseMvc2=_interopRequireDefault(_baseMvc),_localization=require("utils/localization"),_localization2=_interopRequireDefault(_localization),_super=_historyView2.default.HistoryView,AnnotatedHistoryView=_super.extend({className:_super.prototype.className+" annotated-history-panel",_buildNewRender:function(){var e=_super.prototype._buildNewRender.call(this);return this.renderHistoryAnnotation(e),e},renderHistoryAnnotation:function(e){var t=this.model.get("annotation");t&&e.find("> .controls .subtitle").text(t)},renderItems:function(e){e=e||this.$el,_super.prototype.renderItems.call(this,e);var t=e.find("> .controls");t.find(".contents-container.headers").remove();$('<div class="contents-container headers"/>').append([$('<div class="history-content header"/>').text((0,_localization2.default)("Dataset")),$('<div class="additional-info header"/>').text((0,_localization2.default)("Annotation"))]).appendTo(t);return self.views},_renderItemView$el:function(e){return $('<div class="contents-container"/>').append([e.render(0).$el,$('<div class="additional-info"/>').text(e.model.get("annotation")||"")])},events:_.extend(_.clone(_super.prototype.events),{"click .contents-container":function(e){e.stopPropagation(),$(e.currentTarget).find(".list-item .title-bar").click()},"click .icon-btn":function(e){e.stopPropagation();var t=$(e.currentTarget);t.length&&"dropdown"===t.attr("data-toggle")&&t.dropdown("toggle")}}),_clickSectionLink:function(e){var t=$(e.currentTarget).parent().parent().data("section");this.openSection(t)},toString:function(){return"AnnotatedHistoryView("+(this.model?this.model.get("name"):"")+")"}});exports.default={AnnotatedHistoryView:AnnotatedHistoryView};
//# sourceMappingURL=../../../maps/mvc/history/history-view-annotated.js.map
