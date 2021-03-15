let counter = 0;

function count() {
    counter++;
    //alert(counter)
    document.querySelector('#t1').innerHTML = counter;

    if (counter % 10 === 0) {
        alert(`Count is now ${counter}`);
    }
}

function predict() {
    text_li = document.querySelector('li').innerHTML;
    //alert(`${text_list}`);
    if (text_li === 'No tasks.') {
        alert(`Fill in the cell type and rsID first.`);
    } else {

    }

}

document.addEventListener('DOMContentLoaded', function() {
    document.querySelector('#Calculate').onclick = predict; //count;
});