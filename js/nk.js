$(document).ready(function () {


    $('#ParameterInput').on('submit', function (e) {
        e.preventDefault();
        var periods = parseInt($('#periods').val());
        var kappa = parseFloat($('#kappa').val());
        var rho_g = parseFloat($('#rhoG').val());
        var rho_u = parseFloat($('#rhoU').val());
        var sig_g = parseFloat($('#sigEpsG').val());
        var sig_u = parseFloat($('#sigEpsU').val());
        var pi_T = parseFloat($('#piT').val());
        var phi_pi = parseFloat($('#phiPi').val());
        var phi_y = parseFloat($('#phiY').val());


        var beta = 0.995
        var rho_v = 0.5
        var r_bar = 100*(Math.pow(1/beta,4)-1)
        var y_bar = 1
        var sigma = 1



        var a1 = (1-beta*rho_g)/((1-beta*rho_g)*(1-rho_g+phi_y/sigma)+kappa/sigma*(phi_pi-rho_g))
        var a2 = -(phi_pi-rho_u)/sigma/((1-beta*rho_u)*(1-rho_u+phi_y/sigma)+kappa/sigma*(phi_pi-rho_u))
        var a3 = -(1-beta*rho_v)/sigma/((1-beta*rho_v)*(1-rho_v+phi_y/sigma)+kappa/sigma*(phi_pi-rho_v))
        
        var b1 = kappa/((1-beta*rho_g)*(1-rho_g+phi_y/sigma)+kappa/sigma*(phi_pi-rho_g))
        var b2 = (1-rho_u+phi_y/sigma)/((1-beta*rho_u)*(1-rho_u+phi_y/sigma)+kappa/sigma*(phi_pi-rho_u))
        var b3 = -kappa/sigma/((1-beta*rho_v)*(1-rho_v+phi_y/sigma)+kappa/sigma*(phi_pi-rho_v))

        var c1 = phi_y*a1+phi_pi*b1
        var c2 = phi_y*a2+phi_pi*b2
        var c3 = phi_y*a3+phi_pi*b3+1

        var d1 = c1-rho_g*b1
        var d2 = c2-rho_u*b2
        var d3 = c3-rho_v*b3

        var interval = 0;
        var enableMarks = true;


        var z = new ziggurat();


        if (periods<=30) {
            interval = 5;
        }
        else if (30<periods&&periods<=100) {
            interval = 20;
            enableMarks = false;
        }
        else if (100<periods&&periods<=500) {
            interval = 100;
            enableMarks = false;
        }
        else {
            interval = 500;
            enableMarks = false;
        };


        var g = 0;
        var u = 0;
        var v = 0;
        var g_series = [0];
        var u_series = [0];
        var v_series = [0];
        var e_g_series = [0];
        var e_u_series = [0];
        var pi_series = [0];
        var y_series = [0];
        var i_series = [0];
        var r_series = [0];

        
        if (document.getElementById ("irDemand").checked == true) {
            var include_g_shock = 1
        } else {
            var include_g_shock = 0
        }

        if (document.getElementById ("irPi").checked == true) {
            var include_u_shock = 1
        } else {
            var include_u_shock = 0
        }

        if (document.getElementById ("irInterest").checked == true) {
            var include_i_shock = 1
        } else {
            var include_i_shock = 0
        }

        

        var input = document.getElementById ("stochSim");
            if (input.checked == true) {
                for (i = 0; i <= periods; i++) {
                    e_g = z.nextGaussian()
                    e_u = z.nextGaussian()

                    g = rho_g*g+sig_g*e_g
                    u = rho_u*u+sig_u*e_u
                    v = 0

                    e_g_series.push(e_g);
                    e_u_series.push(e_u);

                    g_series.push(100*g);
                    u_series.push(100*u);
                    v_series.push(100*v);

                    y_series.push(100*(a1*g+a2*u+a3*v))
                    pi_series.push(400*(b1*g+b2*u+b3*v))
                    i_series.push(400*(c1*g+c2*u+c3*v))
                    r_series.push(400*(d1*g+d2*u+d3*v))


                }
            } else {
                for (i = 0; i <= periods; i++) {
                    if (i<1) {
                        e_g = sig_g*include_g_shock
                        e_u = sig_u*include_u_shock
                        e_v = 0.01/4*include_i_shock


                        g = rho_g*g+e_g
                        u = rho_u*u+e_u
                        v = rho_v*v+e_v

                        // document.write(g,u,v)

                        e_g_series.push(e_g);
                        e_u_series.push(e_u);

                        g_series.push(100*g);
                        u_series.push(100*u);
                        v_series.push(100*v);

                        y_series.push(100*(a1*g+a2*u+a3*v))
                        pi_series.push(400*(b1*g+b2*u+b3*v))
                        i_series.push(400*(c1*g+c2*u+c3*v))
                        r_series.push(400*(d1*g+d2*u+d3*v))
                    }
                    else {
                        e_g = 0
                        e_u = 0
                        e_v = 0

                        g = rho_g*g+e_g
                        u = rho_u*u+e_u
                        v = rho_v*v+e_v

                        e_g_series.push(e_g);
                        e_u_series.push(e_u);

                        g_series.push(100*g);
                        u_series.push(100*u);
                        v_series.push(100*v);

                        y_series.push(100*(a1*g+a2*u+a3*v))
                        pi_series.push(400*(b1*g+b2*u+b3*v))
                        i_series.push(400*(c1*g+c2*u+c3*v))
                        r_series.push(400*(d1*g+d2*u+d3*v))
                    }

                }
            };

        sessionStorage.setItem('g_series',JSON.stringify(g_series));
        sessionStorage.setItem('u_series',JSON.stringify(u_series));
        sessionStorage.setItem('v_series',JSON.stringify(v_series));
        sessionStorage.setItem('y_series',JSON.stringify(y_series));
        sessionStorage.setItem('pi_series',JSON.stringify(pi_series));
        sessionStorage.setItem('i_series',JSON.stringify(i_series));
        sessionStorage.setItem('r_series',JSON.stringify(r_series));
        sessionStorage.setItem('periods',periods);


        $('#shockProcesses').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'Exogenous Processes (%)',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 0.0001,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: true
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'g',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*g_series[i])/100000,
                            // y: g_series[i],
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'u',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*u_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });

        $('#iProcesses').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'Interest Rates (%)',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 0.0001,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: true
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'i',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*i_series[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            {
                name: 'r',
                color: '#f15c80',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*r_series[i])/100000,
                        })
                    }
                    return data;
                })()
            }
            ]
        });

        $('#piProcess').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'Inflation (%)',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 0.0001,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: false
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'pi',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*pi_series[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            ]
        });
        $('#yProcess').highcharts({
            chart: {
                type: 'line',
                marginRight: 10,
            },

            credits: {
                enabled: false
            },

            plotOptions: {
                series: {
                    lineWidth: 2,
                    marker : {
                        enabled : false,
                        radius : 3},
                    animation: {
                        duration: 10000     //Controls the time required for plot to be fully rendered.
                    }
                }
            },

            title: {
                text: 'Output Gap (%)',
                style: {
                    "fontSize": "15px" 
                }
            },
            xAxis: {
                // type: 'datetime',
                // tickPixelInterval: 150,
                tickInterval: interval,
                text: 'Time'
            },
            yAxis: {
                title: {
                    text: ''
                },
                labels: {
                    formatter: function() 
                    {
                        return Math.round(this.value*100)/100;
                    }
                },
                minRange: 0.0001,
                plotLines: [{
                    value: 0,
                    width: 2,
                    color: '#808080'
                }]
            },

            legend: {
                enabled: false
            },
            exporting: {
                enabled: true,
                buttons: {
                    contextButton: {
                        menuItems: ['downloadPNG']
                    },
                },
            },
            series: [{
                name: 'y',
                data: (function () {
                    var data = [];

                    for (i = 0; i <= periods; i++) {

                        data.push({
                            x: i,
                            y: Math.round(100000*y_series[i])/100000,
                        })
                    }
                    return data;
                })()
            },
            ]
        });
    }) //.trigger('submit');
});

function reloadFunction() {
    sessionStorage.clear();
    location.reload();
}

function downloadFunction() {
    
    if (sessionStorage.getItem('y_series') == null){
        window.alert('Run the simulation first.')
    } else {
        var g_data = JSON.parse(sessionStorage.getItem('g_series'));
        var u_data = JSON.parse(sessionStorage.getItem('u_series'));
        var v_data = JSON.parse(sessionStorage.getItem('v_series'));
        var y_data = JSON.parse(sessionStorage.getItem('y_series'));
        var pi_data = JSON.parse(sessionStorage.getItem('pi_series'));
        var i_data = JSON.parse(sessionStorage.getItem('i_series'));
        var r_data = JSON.parse(sessionStorage.getItem('r_series'));
        var periods = parseInt(sessionStorage.getItem('periods'));

        var row1 = [];
        row1.push("Period");
        row1.push("g");
        row1.push("u");
        row1.push("v");
        row1.push("y");
        row1.push("pi");
        row1.push("i");
        row1.push("r");

        var rows = [row1]
        for (i = 0; i <= periods; i++) {
            rows.push([i,g_data[i],u_data[i],v_data[i],y_data[i],pi_data[i],i_data[i],r_data[i]]);
        }

        let csvContent = "data:text/csv;charset=utf-8,";
        rows.forEach(function(rowArray){
           let row = rowArray.join(",");
           csvContent += row + "\r\n"; // add carriage return
        }); 

        var encodedUri = encodeURI(csvContent);
        var link = document.createElement("a");
        link.setAttribute("href", encodedUri);
        link.setAttribute("download", "nk_data.csv");
        document.body.appendChild(link); // Required for FF

        link.click(); // This will download the data file named "nk_data.csv".
    }

}