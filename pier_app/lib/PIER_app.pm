package PIER_app;
use Mojo::Base 'Mojolicious';

# This method will run once at server start
sub startup {
	my $self = shift;
	
	$ENV{MOJO_REVERSE_PROXY} = 1;
	$self->config(
		hypnotoad => {
			listen  => ['http://*:80'],
			workers => 8,
			keep_alive_timeout => 300,
			websocket_timeout => 600,
			proxy => 1
		}
	);
	
	# Documentation browser under "/perldoc"
  	$self->plugin('PODRenderer');

  	# Router
  	my $r = $self->routes;
	
	# Template names are expected to follow the template.format.handler scheme, with template defaulting to controller/action or the route name, format defaulting to html and handler to ep
	
  	# Normal route to controller
  	## Home
  	$r->get('/')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/PiER')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/pier')->to(template=>'index', controller=>'action', action=>'index');
	$r->get('/PIER')->to(template=>'index', controller=>'action', action=>'index');
	
	##############################
	## booklet (see Action.pm -> sub booklet -> redirect_to)
	## so that '/PiER/booklet' equivalent to '/PiERbooklet/index.html' (located at pier_app/public/PiERbooklet/index.html)
  	$r->get('/PiER/booklet')->to(controller=>'action', action=>'booklet');
	##############################
	
	##############################
  	## PiER
  	### eV2CG
  	$r->get('/PiER/eV2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>0);
  	$r->post('/PiER/eV2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>1);
  	### cTCrosstalk
  	$r->get('/PiER/cTCrosstalk')->to(template=>'cTCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTCrosstalk', post_flag=>0);
  	$r->post('/PiER/cTCrosstalk')->to(template=>'cTCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTCrosstalk', post_flag=>1);
  	### cTGene
  	$r->get('/PiER/cTGene')->to(template=>'cTGene', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTGene', post_flag=>0);
  	$r->post('/PiER/cTGene')->to(template=>'cTGene', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTGene', post_flag=>1);
  	### eV2CG
  	$r->get('/PiER/eV2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>0);
  	$r->post('/PiER/eV2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>1);
  	### eCG2PG
  	$r->get('/PiER/eCG2PG')->to(template=>'eCG2PG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCG2PG', post_flag=>0);
  	$r->post('/PiER/eCG2PG')->to(template=>'eCG2PG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCG2PG', post_flag=>1);
  	### eCrosstalk
  	$r->get('/PiER/eCrosstalk')->to(template=>'eCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCrosstalk', post_flag=>0);
  	$r->post('/PiER/eCrosstalk')->to(template=>'eCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCrosstalk', post_flag=>1);
	
	
	##############################
  	## PiER (old)
  	### eV2CG
  	$r->get('/e/V2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>0);
  	$r->post('/e/V2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>1);
  	### cTCrosstalk
  	$r->get('/c/TCrosstalk')->to(template=>'cTCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTCrosstalk', post_flag=>0);
  	$r->post('/c/TCrosstalk')->to(template=>'cTCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTCrosstalk', post_flag=>1);
  	### cTGene
  	$r->get('/c/TGene')->to(template=>'cTGene', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTGene', post_flag=>0);
  	$r->post('/c/TGene')->to(template=>'cTGene', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_cTGene', post_flag=>1);
  	### eV2CG
  	$r->get('/e/V2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>0);
  	$r->post('/e/V2CG')->to(template=>'eV2CG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eV2CG', post_flag=>1);
  	### eCG2PG
  	$r->get('/e/CG2PG')->to(template=>'eCG2PG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCG2PG', post_flag=>0);
  	$r->post('/e/CG2PG')->to(template=>'eCG2PG', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCG2PG', post_flag=>1);
  	### eCrosstalk
  	$r->get('/e/Crosstalk')->to(template=>'eCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCrosstalk', post_flag=>0);
  	$r->post('/e/Crosstalk')->to(template=>'eCrosstalk', format=>'html', handler=>'ep', controller=>'action', action=>'PiER_eCrosstalk', post_flag=>1);

	##############################
	## FDB
	$r->get('/FDB')->to(controller=>'action', action=>'fdb');
	$r->get('/fdb')->to(controller=>'action', action=>'fdb');
	##############################

}

1;
