"use strict";function _interopRequireDefault(t){return t&&t.__esModule?t:{default:t}}Object.defineProperty(exports,"__esModule",{value:!0});var _graph=require("utils/graph"),_graph2=_interopRequireDefault(_graph),_addLogging=require("utils/add-logging"),_addLogging2=_interopRequireDefault(_addLogging),_super=_graph2.default.Graph,JobDAG=function(t){t=t||{};var o=this;o.filters=[],o._jobsData=[],o._historyContentsMap={},o._toolMap={},o._outputIdToJobMap={},o.noInputJobs=[],o.noOutputJobs=[],o.filteredSetMetadata=[],o.filteredErroredJobs=[],o.dataKeys=["jobs","historyContents","tools"],_super.call(o,!0,_.pick(t,o.dataKeys),_.omit(t,o.dataKeys))};JobDAG.prototype=new _graph2.default.Graph,JobDAG.prototype.constructor=JobDAG,(0,_addLogging2.default)(JobDAG),JobDAG.prototype.init=function(t){t=t||{};var o=this;return o.options=_.defaults(t,{excludeSetMetadata:!1}),o.filters=o._initFilters(),_super.prototype.init.call(o,t),o},JobDAG.prototype._initFilters=function(){var t=this,o=[];return t.options.excludeSetMetadata&&(t.filteredSetMetadata=[],o.push(function(o){return"__SET_METADATA__"!==o.job.tool_id||(t.filteredSetMetadata.push(o.job.id),!1)})),t.options.excludeErroredJobs&&(t.filteredErroredJobs=[],o.push(function(o){return"error"!==o.job.state||(t.filteredErroredJobs.push(o.job.id),!1)})),_.isArray(t.options.filters)&&(o=o.concat(t.options.filters)),t.debug("filters len:",o.length),o},JobDAG.prototype.read=function(t){var o=this;return _.has(t,"historyContents")&&_.has(t,"jobs")&&_.has(t,"tools")?(o.preprocessHistoryContents(t.historyContents||[]).preprocessTools(t.tools||{}).preprocessJobs(t.jobs||[]),o.createGraph(o._filterJobs()),o):_super.prototype.read.call(this,t)},JobDAG.prototype.preprocessHistoryContents=function(t){this.info("processing history");var o=this;return o._historyContentsMap={},t.forEach(function(t,e){o._historyContentsMap[t.id]=_.clone(t)}),o},JobDAG.prototype.preprocessTools=function(t){this.info("processing tools");var o=this;return o._toolMap={},_.each(t,function(t,e){o._toolMap[e]=_.clone(t)}),o},JobDAG.prototype.preprocessJobs=function(t){this.info("processing jobs");var o=this;return o._outputIdToJobMap={},o._jobsData=o.sort(t).map(function(t){return o.preprocessJob(_.clone(t))}),o},JobDAG.prototype.sort=function(t){return t.sort(function(t,o){return t.create_time>o.create_time?1:t.create_time<o.create_time?-1:0})},JobDAG.prototype.preprocessJob=function(t,o){var e=this,r={job:t};return r.inputs=e._processInputs(t),0===_.size(r.inputs)&&e.noInputJobs.push(t.id),r.outputs=e._processOutputs(t),0===_.size(r.outputs)&&e.noOutputJobs.push(t.id),r.tool=e._toolMap[t.tool_id],r},JobDAG.prototype._processInputs=function(t){var o=this,e=t.inputs,r={};return _.each(e,function(t,e){(t=_.clone(o._validateInputOutput(t))).name=e,t.content=o._historyContentsMap[t.id],r[t.id]=t}),r},JobDAG.prototype._validateInputOutput=function(t){if(!t.id)throw new Error("No id on job input/output: ",JSON.stringify(t));if(!t.src||"hda"!==t.src)throw new Error("Bad src on job input/output: ",JSON.stringify(t));return t},JobDAG.prototype._processOutputs=function(t){var o=this,e=t.outputs,r={};return _.each(e,function(e,n){(e=_.clone(o._validateInputOutput(e))).name=n,e.content=o._historyContentsMap[e.id],r[e.id]=e,o._outputIdToJobMap[e.id]=t.id}),r},JobDAG.prototype._filterJobs=function(){var t=this;return t._jobsData.filter(function(o,e){return t._filterJob(o,e)})},JobDAG.prototype._filterJob=function(t,o){for(var e=this,r=0;r<e.filters.length;r++)if(!e.filters[r].call(e,t))return e.debug("\t job",t.job.id," has been filtered out by function:\n",e.filters[r]),!1;return!0},JobDAG.prototype.createGraph=function(t){var o=this;return o.debug("connections:"),_.each(t,function(t){var e=t.job.id;o.debug("\t",e,t),o.createVertex(e,t)}),_.each(t,function(t){var e=t.job.id;_.each(t.inputs,function(t,r){var n=o._outputIdToJobMap[r];n||(n=o.createJobLessVertex(r).name),o.createEdge(n,e,o.directed,{dataset:r})})}),o.debug("final graph: ",JSON.stringify(o.toVerticesAndEdges(),null,"  ")),o},JobDAG.prototype.createJobLessVertex=function(t){var o="copy-"+t;return this.createVertex(o,this._historyContentsMap[t])},JobDAG.prototype.weakComponentGraphArray=function(){var t=this;return this.weakComponents().map(function(o){return o.vertices.sort(function(t,o){var e=t.data.job?t.data.job.create_time:t.data.create_time,r=o.data.job?o.data.job.create_time:o.data.create_time;return e>r?1:e<r?-1:0}),new Graph(t.directed,o)})},JobDAG.prototype._jobsDataMap=function(){var t={};return this._jobsData.forEach(function(o){t[o.job.id]=o}),t},exports.default=JobDAG;
//# sourceMappingURL=../../../maps/mvc/history/job-dag.js.map
