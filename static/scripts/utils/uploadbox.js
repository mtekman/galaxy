"use strict";!function(e){jQuery.event.props.push("dataTransfer"),e.uploadpost=function(n){var r=e.extend({},{data:{},success:function(){},error:function(){},progress:function(){},url:null,maxfilesize:2048,error_filesize:"File exceeds 2GB. Please use a FTP client.",error_default:"Please make sure the file is available.",error_server:"Upload request failed.",error_login:"Uploads require you to log in."},n),t=r.data;if(t.error_message)r.error(t.error_message);else{var o=new FormData;for(var a in t.payload)o.append(a,t.payload[a]);var i=0;for(var a in t.files){var s=t.files[a];o.append(s.name,s.file,s.file.name),i+=s.file.size}if(i>1048576*r.maxfilesize)r.error(r.error_filesize);else{var u=new XMLHttpRequest;u.open("POST",r.url,!0),u.setRequestHeader("Accept","application/json"),u.setRequestHeader("Cache-Control","no-cache"),u.setRequestHeader("X-Requested-With","XMLHttpRequest"),u.onreadystatechange=function(){if(u.readyState==u.DONE){var e=null,n="";if(u.responseText)try{e=jQuery.parseJSON(u.responseText),n=e.err_msg}catch(r){e=u.responseText,n=e}if(u.status<200||u.status>299){var t=u.statusText;403==u.status?t=r.error_login:0==u.status?t=r.error_server:t||(t=r.error_default),r.error(t+" ("+u.status+"). "+n)}else r.success(e)}},u.upload.addEventListener("progress",function(e){e.lengthComputable&&r.progress(Math.round(100*e.loaded/e.total))},!1),Galaxy.emit.debug("uploadbox::uploadpost()","Posting following data.",r),u.send(o)}}},e.fn.uploadinput=function(n){var r=this,t=e.extend({},{ondragover:function(){},ondragleave:function(){},onchange:function(){},multiple:!1},n),o=e('<input type="file" style="display: none" '+(t.multiple&&"multiple"||"")+"/>");return r.append(o.change(function(n){t.onchange(n.target.files),e(this).val("")})),r.on("drop",function(e){t.ondragleave(e),e.dataTransfer&&(t.onchange(e.dataTransfer.files),e.preventDefault())}),r.on("dragover",function(e){e.preventDefault(),t.ondragover(e)}),r.on("dragleave",function(e){e.stopPropagation(),t.ondragleave(e)}),{dialog:function(){o.trigger("click")}}},e.fn.uploadbox=function(n){function r(e){if(e&&e.length&&!l){var n=void 0;return _.each(e,function(e,n){"new"!==e.mode&&_.filter(i,function(n){return n.name===e.name&&n.size===e.size}).length&&(e.duplicate=!0)}),_.each(e,function(e){e.duplicate||(n=String(s++),i[n]=e,a.announce(n,i[n]),u++)}),n}}function t(e){i[e]&&(delete i[e],u--)}function o(){if(0==u||c)return c=!1,l=!1,void a.complete();l=!0;var n=-1;for(var r in i){n=r;break}i[n];t(n),e.uploadpost({url:a.url,data:a.initialize(n),success:function(e){a.success(n,e),o()},error:function(e){a.error(n,e),o()},progress:function(e){a.progress(n,e)}})}var a=e.extend({},{dragover:function(){},dragleave:function(){},announce:function(e){},initialize:function(e){},progress:function(e,n){},success:function(e,n){},error:function(e,n){alert(n)},complete:function(){}},n),i={},s=0,u=0,l=!1,c=!1,f=e(this).uploadinput({multiple:!0,onchange:function(e){r(e)},ondragover:n.ondragover,ondragleave:n.ondragleave});return{select:function(){f.dialog()},add:r,remove:t,start:function(){l||(l=!0,o())},stop:function(){c=!0},reset:function(e){for(e in i)t(e)},configure:function(n){return a=e.extend({},a,n)},compatible:function(){return window.File&&window.FormData&&window.XMLHttpRequest&&window.FileList}}}}(jQuery);
//# sourceMappingURL=../../maps/utils/uploadbox.js.map
