$(document).ready(function () {


    $('#ParameterInput').on('submit', function (e) {
        
        e.preventDefault();
        

        // var a1 = (1-beta*rho_g)/((1-beta*rho_g)*(1-rho_g+phi_y/sigma)+kappa/sigma*(phi_pi-rho_g))
        // var a2 = -(phi_pi-rho_u)/sigma/((1-beta*rho_u)*(1-rho_u+phi_y/sigma)+kappa/sigma*(phi_pi-rho_u))
        

        
        var alpha  = .35
        var beta  = 0.99
        var delta   = 0.025
        var rhoa = .9
        var sigma = 1.5
        var phi = 1.5
        var A = 1

        var klratio = Math.pow((A*alpha/(Math.pow(beta,-1)+delta-1)),1/(1-alpha))

        
        // Optimization example. First, define an objective function.

        var rhs = Math.pow(A*klratio**alpha-delta*klratio,-sigma)*(1-alpha)*A*Math.pow(klratio,alpha)

        // document.write(rhs)

        sqr = function(x) { return x*x; };
        objective = function(x) { return Math.pow(rhs - phi/(1-x)*Math.pow(x,sigma),2); }
        // objective = function(x) { rhs ; }

        ff = objective(0.37071457143060327)
        result = numeric.uncmin(objective,[.5])

        l_steady_state = result['solution'][0]

        // 0.37071457143060327

        // document.write(rhs-ff)
        // document.write(ff["solution"])
        console.log(l_steady_state)
        
        // sessionStorage.setItem('r_series',JSON.stringify(r_series));
        // sessionStorage.setItem('periods',periods);


        
    }) //.trigger('submit');
});

function reloadFunction() {
    sessionStorage.clear();
    location.reload();
}

// function downloadFunction() {
    
//     if (sessionStorage.getItem('y_series') == null){
//         window.alert('Run the simulation first.')
//     } else {
//         var g_data = JSON.parse(sessionStorage.getItem('g_series'));
//         var u_data = JSON.parse(sessionStorage.getItem('u_series'));
//         var v_data = JSON.parse(sessionStorage.getItem('v_series'));
//         var y_data = JSON.parse(sessionStorage.getItem('y_series'));
//         var pi_data = JSON.parse(sessionStorage.getItem('pi_series'));
//         var i_data = JSON.parse(sessionStorage.getItem('i_series'));
//         var r_data = JSON.parse(sessionStorage.getItem('r_series'));
//         var periods = parseInt(sessionStorage.getItem('periods'));

//         var row1 = [];
//         row1.push("Period");
//         row1.push("g");
//         row1.push("u");
//         row1.push("v");
//         row1.push("y");
//         row1.push("pi");
//         row1.push("i");
//         row1.push("r");

//         var rows = [row1]
//         for (i = 0; i <= periods; i++) {
//             rows.push([i,g_data[i],u_data[i],v_data[i],y_data[i],pi_data[i],i_data[i],r_data[i]]);
//         }

//         let csvContent = "data:text/csv;charset=utf-8,";
//         rows.forEach(function(rowArray){
//            let row = rowArray.join(",");
//            csvContent += row + "\r\n"; // add carriage return
//         }); 

//         var encodedUri = encodeURI(csvContent);
//         var link = document.createElement("a");
//         link.setAttribute("href", encodedUri);
//         link.setAttribute("download", "nk_data.csv");
//         document.body.appendChild(link); // Required for FF

//         link.click(); // This will download the data file named "nk_data.csv".
//     }

// }