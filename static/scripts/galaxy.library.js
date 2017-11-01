"use strict";define(["layout/masthead","utils/utils","libs/toastr","mvc/base-mvc","mvc/library/library-model","mvc/library/library-folderlist-view","mvc/library/library-librarylist-view","mvc/library/library-librarytoolbar-view","mvc/library/library-foldertoolbar-view","mvc/library/library-dataset-view","mvc/library/library-library-view","mvc/library/library-folder-view"],function(r,i,e,a,l,o,s,t,n,d,b,w){var y=Backbone.Router.extend({initialize:function(){this.routesHit=0,Backbone.history.on("route",function(){this.routesHit++},this),this.bind("route",this.trackPageview)},routes:{"":"libraries","page/:show_page":"libraries_page","library/:library_id/permissions":"library_permissions","folders/:folder_id/permissions":"folder_permissions","folders/:id":"folder_content","folders/:id/page/:show_page":"folder_page","folders/:folder_id/datasets/:dataset_id":"dataset_detail","folders/:folder_id/datasets/:dataset_id/permissions":"dataset_permissions","folders/:folder_id/datasets/:dataset_id/versions/:ldda_id":"dataset_version","folders/:folder_id/download/:format":"download","folders/:folder_id/import/:source":"import_datasets"},back:function(){this.routesHit>1?window.history.back():this.navigate("#",{trigger:!0,replace:!0})},trackPageview:function(){var r=Backbone.history.getFragment();/^\//.test(r)||""==r||(r="/"+r),"undefined"!=typeof ga&&ga("send","pageview",Galaxy.root+"library/list"+r)}}),u=a.SessionStorageModel.extend({defaults:{with_deleted:!1,without_restricted:!1,sort_order:"asc",sort_by:"name",library_page_size:20,folder_page_size:15}});return{GalaxyApp:Backbone.View.extend({libraryToolbarView:null,libraryListView:null,library_router:null,libraryView:null,folderToolbarView:null,folderListView:null,datasetView:null,initialize:function(){window.Galaxy.config.ga_code&&(!function(r,i,e,a,l,o,s){r.GoogleAnalyticsObject=l,r[l]=r[l]||function(){(r[l].q=r[l].q||[]).push(arguments)},r[l].l=1*new Date,o=i.createElement(e),s=i.getElementsByTagName(e)[0],o.async=1,o.src="//www.google-analytics.com/analytics.js",s.parentNode.insertBefore(o,s)}(window,document,"script",0,"ga"),ga("create",window.Galaxy.config.ga_code,"auto"),ga("send","pageview")),Galaxy.libraries=this,this.preferences=new u({id:"global-lib-prefs"}),this.library_router=new y,this.library_router.on("route:libraries",function(){Galaxy.libraries.libraryToolbarView&&Galaxy.libraries.libraryToolbarView.$el.unbind("click"),Galaxy.libraries.libraryToolbarView=new t.LibraryToolbarView,Galaxy.libraries.libraryListView=new s.LibraryListView}),this.library_router.on("route:libraries_page",function(r){null===Galaxy.libraries.libraryToolbarView?(Galaxy.libraries.libraryToolbarView=new t.LibraryToolbarView,Galaxy.libraries.libraryListView=new s.LibraryListView({show_page:r})):Galaxy.libraries.libraryListView.render({show_page:r})}),this.library_router.on("route:folder_content",function(r){Galaxy.libraries.folderToolbarView&&Galaxy.libraries.folderToolbarView.$el.unbind("click"),Galaxy.libraries.folderToolbarView=new n.FolderToolbarView({id:r}),Galaxy.libraries.folderListView=new o.FolderListView({id:r})}),this.library_router.on("route:folder_page",function(r,i){null===Galaxy.libraries.folderToolbarView?(Galaxy.libraries.folderToolbarView=new n.FolderToolbarView({id:r}),Galaxy.libraries.folderListView=new o.FolderListView({id:r,show_page:i})):Galaxy.libraries.folderListView.render({id:r,show_page:parseInt(i)})}),this.library_router.on("route:download",function(r,i){0===$("#folder_list_body").find(":checked").length?(e.info("You must select at least one dataset to download"),Galaxy.libraries.library_router.navigate("folders/"+r,{trigger:!0,replace:!0})):(Galaxy.libraries.folderToolbarView.download(r,i),Galaxy.libraries.library_router.navigate("folders/"+r,{trigger:!1,replace:!0}))}),this.library_router.on("route:dataset_detail",function(r,i){Galaxy.libraries.datasetView&&Galaxy.libraries.datasetView.$el.unbind("click"),Galaxy.libraries.datasetView=new d.LibraryDatasetView({id:i,show_version:!1,show_permissions:!1})}),this.library_router.on("route:dataset_version",function(r,i,e){Galaxy.libraries.datasetView&&Galaxy.libraries.datasetView.$el.unbind("click"),Galaxy.libraries.datasetView=new d.LibraryDatasetView({id:i,ldda_id:e,show_version:!0})}),this.library_router.on("route:dataset_permissions",function(r,i){Galaxy.libraries.datasetView&&Galaxy.libraries.datasetView.$el.unbind("click"),Galaxy.libraries.datasetView=new d.LibraryDatasetView({id:i,show_permissions:!0})}),this.library_router.on("route:library_permissions",function(r){Galaxy.libraries.libraryView&&Galaxy.libraries.libraryView.$el.unbind("click"),Galaxy.libraries.libraryView=new b.LibraryView({id:r,show_permissions:!0})}),this.library_router.on("route:folder_permissions",function(r){Galaxy.libraries.folderView&&Galaxy.libraries.folderView.$el.unbind("click"),Galaxy.libraries.folderView=new w.FolderView({id:r,show_permissions:!0})}),this.library_router.on("route:import_datasets",function(r,i){Galaxy.libraries.folderToolbarView&&Galaxy.libraries.folderListView?Galaxy.libraries.folderToolbarView.showImportModal({source:i}):(Galaxy.libraries.folderToolbarView=new n.FolderToolbarView({id:r}),Galaxy.libraries.folderListView=new o.FolderListView({id:r}),Galaxy.libraries.folderToolbarView.showImportModal({source:i}))}),Backbone.history.start({pushState:!1})}})}});
//# sourceMappingURL=../maps/galaxy.library.js.map
